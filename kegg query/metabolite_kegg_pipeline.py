#!/usr/bin/env python3
"""
Metabolomics name -> KEGG candidate mapping pipeline

Inputs:
  - CSV with a metabolite name column (default: "name")

Outputs:
  1) normalized.jsonl: one JSON record per input name (LLM-like fields)
  2) kegg_candidates.jsonl: KEGG candidates with scores + decision
  3) ranked.csv: flat table (one row per candidate) + per-name decision

Notes:
  - This script uses deterministic heuristics for "LLM-like" 
normalization.
  - It uses KEGG REST for candidate generation and DOES NOT invent KEGG 
IDs.
  - Rate-limits KEGG calls and caches responses locally to avoid hammering 
KEGG.
"""

import argparse
import csv
import json
import os
import re
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote

import requests

try:
    from rapidfuzz import fuzz
    def similarity(a: str, b: str) -> float:
        return float(fuzz.token_set_ratio(a, b)) / 100.0
except Exception:
    # fallback if rapidfuzz isn't installed
    import difflib
    def similarity(a: str, b: str) -> float:
        return difflib.SequenceMatcher(None, a, b).ratio()


# --------------------------
# 1) "LLM-like" JSON schema
# --------------------------
@dataclass
class NormalizedRecord:
    input: str
    normalized: str
    flags: List[str]
    synonyms: List[str]
    must_disambiguate: List[Dict[str, Any]]


# --------------------------
# String normalization
# --------------------------
_STEREO_PREFIX = re.compile(r"^(?:[dDlL]|[rRsS]|[αβ]|alpha|beta)\s*[- ]\s*", re.IGNORECASE)
_MULTI_WS = re.compile(r"\s+")
_BRACKETS = re.compile(r"[\[\]{}()]")
_PUNCT = re.compile(r"[,;]+")

# hydrate / salt-ish tokens that often break text searches
_HYDRATE = re.compile(r"\b(monohydrate|dihydrate|trihydrate|tetrahydrate|hydrate)\b", re.IGNORECASE)

# crude positional-isomer ambiguity triggers
# e.g. phenylbutyrate (but not 2-phenylbutyrate / 4-phenylbutyrate)
_POS_ISOMER_PATTERNS = [
    (re.compile(r"\bphenylbutyrate\b", re.IGNORECASE), 
["2-phenylbutyrate", "3-phenylbutyrate", "4-phenylbutyrate"]),
    (re.compile(r"\bphenylbutyric acid\b", re.IGNORECASE), 
["2-phenylbutyric acid", "3-phenylbutyric acid", "4-phenylbutyric acid"]),
]

_POS_SPECIFIED = re.compile(r"\b[2-6]\s*-\s*", re.IGNORECASE)


def normalize_name(raw: str) -> str:
    s = raw.strip()
    s = _PUNCT.sub(" ", s)
    s = _BRACKETS.sub(" ", s)
    s = _HYDRATE.sub("", s)
    s = _STEREO_PREFIX.sub("", s)  # drop leading stereo token (keeps a "looser" normalized form)
    s = _MULTI_WS.sub(" ", s).strip()
    return s


def synonym_expansion(raw: str, normalized: str) -> List[str]:
    syns = []
    r = raw.strip()

    # basic variants
    if normalized.lower() != r.lower():
        syns.append(r)

    # remove/replace common separators
    syns.append(normalized.replace("-", " "))
    syns.append(normalized.replace(" ", "-"))

    # drop commas already handled; also add lowercase variant
    syns.append(normalized.lower())

    # if ends with "acid" consider "ate" (very heuristic)
    if normalized.lower().endswith(" acid"):
        syns.append(re.sub(r"\s+acid$", "ate", normalized, flags=re.IGNORECASE))

    # remove duplicate/empty
    out = []
    seen = set()
    for x in syns:
        x2 = _MULTI_WS.sub(" ", x).strip()
        if not x2:
            continue
        key = x2.lower()
        if key not in seen:
            seen.add(key)
            out.append(x2)
    return out[:25]


def build_normalized_record(name: str) -> NormalizedRecord:
    norm = normalize_name(name)
    flags: List[str] = []
    must: List[Dict[str, Any]] = []

    # positional isomer ambiguity (only if not already specified like "2-...")
    if not _POS_SPECIFIED.search(name):
        for pat, candidates in _POS_ISOMER_PATTERNS:
            if pat.search(name):
                flags.append("positional_isomer_ambiguous")
                must.append({"type": "positional_isomer", "candidates": candidates})
                break

    # if hydrate in raw, flag it (useful provenance)
    if _HYDRATE.search(name):
        flags.append("contains_hydrate_token_removed")

    syns = synonym_expansion(name, norm)
    return NormalizedRecord(
        input=name,
        normalized=norm,
        flags=flags,
        synonyms=syns,
        must_disambiguate=must
    )


# --------------------------
# 2) KEGG REST querying (cached + rate-limited)
# --------------------------
KEGG_FIND = "https://rest.kegg.jp/find/compound/{}"
KEGG_GET  = "https://rest.kegg.jp/get/{}"

def load_cache(cache_path: Path) -> Dict[str, Any]:
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    return {}

def save_cache(cache_path: Path, cache: Dict[str, Any]) -> None:
    cache_path.write_text(json.dumps(cache))

def kegg_find_compound(query: str, session: requests.Session, cache: 
Dict[str, Any], callsleep: float) -> List[Dict[str, str]]:
    q = query.strip()
    ck = f"find::{q.lower()}"
    if ck in cache:
        return cache[ck]

    url = KEGG_FIND.format(quote(q, safe=""))
    r = session.get(url, timeout=30)
    time.sleep(callsleep)

    hits: List[Dict[str, str]] = []
    if r.status_code == 200 and r.text.strip():
        for line in r.text.strip().splitlines():
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            kegg_id, desc = parts
            hits.append({"kegg_id": kegg_id, "desc": desc})
    cache[ck] = hits
    return hits

def parse_kegg_get(text: str) -> Dict[str, Any]:
    # KEGG flatfile parse for NAME / FORMULA / EXACT_MASS
    out: Dict[str, Any] = {"names": [], "formula": None, "exact_mass": 
None}
    current_key = None
    for line in text.splitlines():
        if not line.strip():
            continue
        if line[:12].strip():  # new field
            key = line[:12].strip()
            val = line[12:].strip()
            current_key = key
        else:
            key = current_key
            val = line[12:].strip()

        if key == "NAME":
            # NAME field may have multiple entries separated by ';'
            val = val.rstrip(";")
            out["names"].append(val)
        elif key == "FORMULA":
            out["formula"] = val
        elif key == "EXACT_MASS":
            try:
                out["exact_mass"] = float(val)
            except Exception:
                out["exact_mass"] = None
    # de-dupe names
    seen = set()
    uniq = []
    for n in out["names"]:
        k = n.lower()
        if k not in seen:
            seen.add(k)
            uniq.append(n)
    out["names"] = uniq
    return out

def kegg_get_compound(kegg_id: str, session: requests.Session, cache: 
Dict[str, Any], callsleep: float) -> Dict[str, Any]:
    ck = f"get::{kegg_id}"
    if ck in cache:
        return cache[ck]
    url = KEGG_GET.format(quote(kegg_id, safe=":"))
    r = session.get(url, timeout=30)
    time.sleep(callsleep)
    if r.status_code != 200:
        cache[ck] = {}
        return {}
    parsed = parse_kegg_get(r.text)
    cache[ck] = parsed
    return parsed


# --------------------------
# 3) Ranking + decision
# --------------------------
def score_candidate(norm: str, syns: List[str], kegg_desc: str, 
kegg_names: List[str]) -> float:
    # Score against KEGG description and NAME field(s)
    targets = [kegg_desc] + (kegg_names or [])
    best = 0.0
    for t in targets:
        best = max(best, similarity(norm.lower(), t.lower()))
        for s in syns[:10]:
            best = max(best, similarity(s.lower(), t.lower()))
    return best

def decide_from_scores(scores: List[float], top_delta: float = 0.07, 
exact_thresh: float = 0.90) -> str:
    if not scores:
        return "no_hit"
    s1 = scores[0]
    s2 = scores[1] if len(scores) > 1 else 0.0
    if s1 >= exact_thresh and (s1 - s2) >= top_delta:
        return "exact"
    return "ambiguous"

def pipeline(
    csv_path: Path,
    out_prefix: Path,
    name_col: str = "name",
    top_k: int = 10,
    calls_per_sec: float = 3.0,
) -> Tuple[Path, Path, Path]:

    callsleep = max(0.0, 1.0 / calls_per_sec)

    cache_path = out_prefix.with_suffix(".kegg_cache.json")
    cache = load_cache(cache_path)

    normalized_path = out_prefix.with_suffix(".normalized.jsonl")
    candidates_path = out_prefix.with_suffix(".kegg_candidates.jsonl")
    ranked_csv_path = out_prefix.with_suffix(".kegg_ranked.csv")

    # read CSV
    rows = []
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if name_col not in reader.fieldnames:
            raise ValueError(f"Column '{name_col}' not found. Found: {reader.fieldnames}")
        for r in reader:
            nm = (r.get(name_col) or "").strip()
            if nm:
                rows.append(nm)

    # de-dupe, preserve order
    seen = set()
    names = []
    for nm in rows:
        k = nm.lower()
        if k not in seen:
            seen.add(k)
            names.append(nm)

    session = requests.Session()
    session.headers.update({"User-Agent": "metabolite-kegg-linker/0.1 (academic use)"})

    # write outputs
    with normalized_path.open("w", encoding="utf-8") as f_norm, \
         candidates_path.open("w", encoding="utf-8") as f_cand, \
         ranked_csv_path.open("w", newline="", encoding="utf-8") as f_rank:

        rank_writer = csv.DictWriter(
            f_rank,
            fieldnames=[
                "input", "normalized", "decision",
                "kegg_id", "kegg_desc", "kegg_name_best",
                "score", "formula", "exact_mass",
                "flags", "match_note"
            ],
        )
        rank_writer.writeheader()

        for nm in names:
            rec = build_normalized_record(nm)
            f_norm.write(json.dumps(asdict(rec), ensure_ascii=False) + 
"\n")

            # generate KEGG candidates using normalized + synonyms
            queries = [rec.normalized] + [s for s in rec.synonyms if s.lower() != rec.normalized.lower()]
            queries = queries[:10]  # token/call discipline

            raw_hits: Dict[str, Dict[str, str]] = {}
            for q in queries:
                hits = kegg_find_compound(q, session, cache, callsleep)
                for h in hits:
                    raw_hits[h["kegg_id"]] = {"kegg_id": h["kegg_id"], "desc": h["desc"], "from_query": q}

            # enrich with KEGG GET for top hits (optional but useful)
            cand_list = []
            for kid, meta in raw_hits.items():
                info = kegg_get_compound(kid, session, cache, callsleep)
                names_kegg = info.get("names", []) if isinstance(info, dict) else []
                formula = info.get("formula") if isinstance(info, dict) else None
                exact_mass = info.get("exact_mass") if isinstance(info, dict) else None

                sc = score_candidate(rec.normalized, rec.synonyms, meta["desc"], names_kegg)
                best_name = names_kegg[0] if names_kegg else None

                cand_list.append({
                    "db": "KEGG",
                    "id": kid,
                    "name": best_name,
                    "desc": meta["desc"],
                    "from_query": meta["from_query"],
                    "formula": formula,
                    "exact_mass": exact_mass,
                    "score": round(sc, 4),
                })

            cand_list.sort(key=lambda x: x["score"], reverse=True)
            top = cand_list[:top_k]

            scores = [c["score"] for c in top]
            decision = decide_from_scores(scores)

            # override: if must_disambiguate triggered, stay ambiguous unless perfect single hit
            if rec.must_disambiguate and decision == "exact":
                # require very strong evidence to override
                if not (scores and scores[0] >= 0.97):
                    decision = "ambiguous"

            result = {
                "input": rec.input,
                "normalized": rec.normalized,
                "flags": rec.flags,
                "must_disambiguate": rec.must_disambiguate,
                "decision": decision,
                "top_candidates": top,
                "missing_info_to_resolve": (["positional/stereo isomer ambiguity; need MS/MS, RT standard, or explicit name"]
                    if decision == "ambiguous" and rec.must_disambiguate
                    else []
                ),
                "recommended_next_step": (["use MS/MS library match", "use RT against authentic standard", "keep ambiguous set or class-level annotation"]
                    if decision == "ambiguous"
                    else []
                ),
            }
            f_cand.write(json.dumps(result, ensure_ascii=False) + "\n")

            # write flattened rows to CSV
            if not top:
                rank_writer.writerow({
                    "input": rec.input,
                    "normalized": rec.normalized,
                    "decision": "no_hit",
                    "kegg_id": "",
                    "kegg_desc": "",
                    "kegg_name_best": "",
                    "score": "",
                    "formula": "",
                    "exact_mass": "",
                    "flags": ";".join(rec.flags),
                    "match_note": "no candidates returned by KEGG find",
                })
            else:
                for c in top:
                    rank_writer.writerow({
                        "input": rec.input,
                        "normalized": rec.normalized,
                        "decision": decision,
                        "kegg_id": c["id"],
                        "kegg_desc": c["desc"],
                        "kegg_name_best": c["name"] or "",
                        "score": c["score"],
                        "formula": c["formula"] or "",
                        "exact_mass": c["exact_mass"] if c["exact_mass"] is not None else "",
                        "flags": ";".join(rec.flags),
                        "match_note": f"from_query={c['from_query']}",
                    })

    save_cache(cache_path, cache)
    return normalized_path, candidates_path, ranked_csv_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, type=Path, help="Input CSV (e.g., mycpds_short.csv)")
    ap.add_argument("--out", required=True, type=Path, help="Output prefix path (no extension)")
    ap.add_argument("--name-col", default="name", help="Column containing metabolite names (default: name)")
    ap.add_argument("--top-k", type=int, default=10, help="Top K candidates to keep per metabolite")
    ap.add_argument("--calls-per-sec", type=float, default=3.0, help="Rate limit for KEGG calls (~3/sec recommended)")
    args = ap.parse_args()

    npath, cpath, rpath = pipeline(
        csv_path=args.csv,
        out_prefix=args.out,
        name_col=args.name_col,
        top_k=args.top_k,
        calls_per_sec=args.calls_per_sec,
    )
    print(f"Wrote: {npath}")
    print(f"Wrote: {cpath}")
    print(f"Wrote: {rpath}")

if __name__ == "__main__":
    main()

