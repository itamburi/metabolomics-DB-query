#!/usr/bin/env python3
"""
Metabolomics name -> HMDB candid:contentReference[oaicite:4]{index=4}
Modeled after metabolite_kegg_pipeline.py (structure preserved):
  - normalized.jsonl: per-input normalization record
  - hmdb_candidates.jsonl: candidates + decision
  - hmdb_ranked.csv: flat table of candidates

HMDB caveat:
- HMDB states official API access requires contacting them.
  So this script uses public web endpoints:
    (1) Candidate generation: /unearth/q (text search)
    (2) Candidate enrichment: metabolite XML, else HTML fallback

Be polite: use caching + rate limiting.
"""

import argparse
import csv
import json
import re
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote_plus

import requests

try:
    from rapidfuzz import fuzz
    def similarity(a: str, b: str) -> float:
        return float(fuzz.token_set_ratio(a, b)) / 100.0
except Exception:
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
# String normalization (copied from KEGG script)
# --------------------------
_STEREO_PREFIX = re.compile(r"^(?:[dDlL]|[rRsS]|[αβ]|alpha|beta)\s*[- ]\s*", re.IGNORECASE)
_MULTI_WS = re.compile(r"\s+")
_BRACKETS = re.compile(r"[\[\]{}()]")
_PUNCT = re.compile(r"[,;]+")
_HYDRATE = re.compile(r"\b(monohydrate|dihydrate|trihydrate|tetrahydrate|hydrate)\b", re.IGNORECASE)

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
    s = _STEREO_PREFIX.sub("", s)
    s = _MULTI_WS.sub(" ", s).strip()
    return s


def synonym_expansion(raw: str, normalized: str) -> List[str]:
    syns = []
    r = raw.strip()

    if normalized.lower() != r.lower():
        syns.append(r)

    syns.append(normalized.replace("-", " "))
    syns.append(normalized.replace(" ", "-"))
    syns.append(normalized.lower())

    if normalized.lower().endswith(" acid"):
        syns.append(re.sub(r"\s+acid$", "ate", normalized, flags=re.IGNORECASE))

    out, seen = [], set()
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

    if not _POS_SPECIFIED.search(name):
        for pat, candidates in _POS_ISOMER_PATTERNS:
            if pat.search(name):
                flags.append("positional_isomer_ambiguous")
                must.append({"type": "positional_isomer", "candidates": candidates})
                break

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
# 2) HMDB web querying (cached + rate-limited)
# --------------------------
HMDB_UNEARTH = "https://www.hmdb.ca/unearth/q?query={query}&searcher=metabolites"
HMDB_XML = "https://hmdb.ca/metabolites/{hmdb_id}.xml"
HMDB_HTML = "https://hmdb.ca/metabolites/{hmdb_id}"

_HMDB_ID_RE = re.compile(r"\bHMDB\d{7}\b")


def load_cache(cache_path: Path) -> Dict[str, Any]:
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    return {}


def save_cache(cache_path: Path, cache: Dict[str, Any]) -> None:
    cache_path.write_text(json.dumps(cache))


def hmdb_find_compound(query: str, session: requests.Session, cache: Dict[str, Any], callsleep: float) -> List[Dict[str, str]]:
    """
    Analogue of KEGG /find/compound:
    - Hits HMDB /unearth/q search
    - Extracts HMDB IDs from HTML
    """
    q = query.strip()
    ck = f"find::{q.lower()}"
    if ck in cache:
        return cache[ck]

    url = HMDB_UNEARTH.format(query=quote_plus(q))
    r = session.get(url, timeout=30)
    time.sleep(callsleep)

    hits: List[Dict[str, str]] = []
    if r.status_code == 200 and r.text:
        ids = _HMDB_ID_RE.findall(r.text)
        seen = set()
        for hmdb_id in ids:
            if hmdb_id in seen:
                continue
            seen.add(hmdb_id)
            hits.append({"hmdb_id": hmdb_id, "desc": ""})

    cache[ck] = hits
    return hits


def _tag_suffix(el: ET.Element) -> str:
    return el.tag.split("}")[-1] if "}" in el.tag else el.tag


def _first_text_by_tags(root: ET.Element, tags: List[str]) -> Optional[str]:
    for el in root.iter():
        if _tag_suffix(el) in tags and el.text and el.text.strip():
            return el.text.strip()
    return None


def _all_text_by_tag(root: ET.Element, tag: str) -> List[str]:
    out = []
    for el in root.iter():
        if _tag_suffix(el) == tag and el.text and el.text.strip():
            out.append(el.text.strip())
    # de-dupe preserve order
    uniq, seen = [], set()
    for x in out:
        k = x.lower()
        if k not in seen:
            seen.add(k)
            uniq.append(x)
    return uniq


def parse_hmdb_xml(text: str) -> Dict[str, Any]:
    """
    Best-effort HMDB XML parsing. HMDB XML structure can vary by version/entry,
    so we do namespace-agnostic tag suffix matching.
    """
    out: Dict[str, Any] = {"names": [], "formula": None, "exact_mass": None}

    root = ET.fromstring(text)

    # Names
    common = _first_text_by_tags(root, ["name", "common_name", "metabolite_name"])
    syns = _all_text_by_tag(root, "synonym")

    names: List[str] = []
    if common and not common.startswith("HMDB"):
        names.append(common)
    names.extend(syns)

    # Formula
    out["formula"] = _first_text_by_tags(root, ["chemical_formula", "formula"])

    # Mass (HMDB often stores monoisotopic as "monisotopic_molecular_weight")
    mono = _first_text_by_tags(root, ["monisotopic_molecular_weight", "monoisotopic_molecular_weight", "exact_mass"])
    try:
        out["exact_mass"] = float(mono) if mono is not None else None
    except Exception:
        out["exact_mass"] = None

    # De-dupe names
    uniq, seen = [], set()
    for n in names:
        k = n.lower()
        if k not in seen:
            seen.add(k)
            uniq.append(n)
    out["names"] = uniq[:50]
    return out


def hmdb_get_compound(hmdb_id: str, session: requests.Session, cache: Dict[str, Any], callsleep: float) -> Dict[str, Any]:
    """
    Analogue of KEGG /get/{id}:
    - Tries HMDB XML metabocard
    - Falls back to minimal HTML-based name extraction
    """
    ck = f"get::{hmdb_id}"
    if ck in cache:
        return cache[ck]

    # XML first
    url_xml = HMDB_XML.format(hmdb_id=hmdb_id)
    r = session.get(url_xml, timeout=30, headers={"Accept": "application/xml,text/xml,*/*"})
    time.sleep(callsleep)

    if r.status_code == 200 and r.text and "<" in r.text[:100]:
        try:
            parsed = parse_hmdb_xml(r.text)
            cache[ck] = parsed
            return parsed
        except Exception:
            pass

    # HTML fallback (very light)
    url_html = HMDB_HTML.format(hmdb_id=hmdb_id)
    r2 = session.get(url_html, timeout=30)
    time.sleep(callsleep)

    if r2.status_code == 200 and r2.text:
        m = re.search(r"Showing metabocard for\s+([^<\(]+)\s*\(", r2.text, flags=re.IGNORECASE)
        nm = m.group(1).strip() if m else None
        parsed = {"names": [nm] if nm else [], "formula": None, "exact_mass": None}
        cache[ck] = parsed
        return parsed

    cache[ck] = {}
    return {}


# --------------------------
# 3) Ranking + decision (same logic as KEGG script)
# --------------------------
def score_candidate(norm: str, syns: List[str], desc: str, names: List[str]) -> float:
    targets = [desc] if desc else []
    targets += (names or [])
    best = 0.0
    for t in targets:
        best = max(best, similarity(norm.lower(), t.lower()))
        for s in syns[:10]:
            best = max(best, similarity(s.lower(), t.lower()))
    return best


def decide_from_scores(scores: List[float], top_delta: float = 0.07, exact_thresh: float = 0.90) -> str:
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
    calls_per_sec: float = 2.0,
) -> Tuple[Path, Path, Path]:

    callsleep = max(0.0, 1.0 / calls_per_sec)

    cache_path = out_prefix.with_suffix(".hmdb_cache.json")
    cache = load_cache(cache_path)

    normalized_path = out_prefix.with_suffix(".normalized.jsonl")
    candidates_path = out_prefix.with_suffix(".hmdb_candidates.jsonl")
    ranked_csv_path = out_prefix.with_suffix(".hmdb_ranked.csv")

    # read CSV
    rows = []
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if name_col not in (reader.fieldnames or []):
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
    session.headers.update({"User-Agent": "metabolite-hmdb-linker/0.1 (academic use; rate-limited)"})

    with normalized_path.open("w", encoding="utf-8") as f_norm, \
         candidates_path.open("w", encoding="utf-8") as f_cand, \
         ranked_csv_path.open("w", newline="", encoding="utf-8") as f_rank:

        rank_writer = csv.DictWriter(
            f_rank,
            fieldnames=[
                "input", "normalized", "decision",
                "hmdb_id", "hmdb_name_best", "score",
                "formula", "exact_mass",
                "flags", "match_note"
            ],
        )
        rank_writer.writeheader()

        for nm in names:
            rec = build_normalized_record(nm)
            f_norm.write(json.dumps(asdict(rec), ensure_ascii=False) + "\n")

            queries = [rec.normalized] + [s for s in rec.synonyms if s.lower() != rec.normalized.lower()]
            queries = queries[:10]

            raw_hits: Dict[str, Dict[str, str]] = {}
            for q in queries:
                hits = hmdb_find_compound(q, session, cache, callsleep)
                for h in hits:
                    raw_hits[h["hmdb_id"]] = {"hmdb_id": h["hmdb_id"], "desc": h.get("desc", ""), "from_query": q}

            cand_list = []
            for hid, meta in raw_hits.items():
                info = hmdb_get_compound(hid, session, cache, callsleep)
                names_hmdb = info.get("names", []) if isinstance(info, dict) else []
                formula = info.get("formula") if isinstance(info, dict) else None
                exact_mass = info.get("exact_mass") if isinstance(info, dict) else None

                sc = score_candidate(rec.normalized, rec.synonyms, meta["desc"], names_hmdb)
                best_name = names_hmdb[0] if names_hmdb else None

                cand_list.append({
                    "db": "HMDB",
                    "id": hid,
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

            if rec.must_disambiguate and decision == "exact":
                if not (scores and scores[0] >= 0.97):
                    decision = "ambiguous"

            result = {
                "input": rec.input,
                "normalized": rec.normalized,
                "flags": rec.flags,
                "must_disambiguate": rec.must_disambiguate,
                "decision": decision,
                "top_candidates": top,
                "missing_info_to_resolve": ([
                    "positional/stereo isomer ambiguity; need MS/MS, RT standard, or explicit name"
                ] if decision == "ambiguous" and rec.must_disambiguate else []),
                "recommended_next_step": ([
                    "use MS/MS library match",
                    "use RT against authentic standard",
                    "keep ambiguous set or class-level annotation"
                ] if decision == "ambiguous" else []),
            }
            f_cand.write(json.dumps(result, ensure_ascii=False) + "\n")

            if not top:
                rank_writer.writerow({
                    "input": rec.input,
                    "normalized": rec.normalized,
                    "decision": "no_hit",
                    "hmdb_id": "",
                    "hmdb_name_best": "",
                    "score": "",
                    "formula": "",
                    "exact_mass": "",
                    "flags": ";".join(rec.flags),
                    "match_note": "no candidates returned by HMDB unearth search",
                })
            else:
                for c in top:
                    rank_writer.writerow({
                        "input": rec.input,
                        "normalized": rec.normalized,
                        "decision": decision,
                        "hmdb_id": c["id"],
                        "hmdb_name_best": c["name"] or "",
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
    ap.add_argument("--calls-per-sec", type=float, default=2.0, help="Rate limit for HMDB calls (be polite)")
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

