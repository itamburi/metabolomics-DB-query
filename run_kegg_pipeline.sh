#!/usr/bin/env bash
set -euo pipefail

# Activate conda env
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate metabolite-query

# Run pipeline
python metabolite_kegg_pipeline.py \
  --csv mycpds_short.csv \
  --out mycpds_short_out \
  --name-col name \
  --top-k 10 \
  --calls-per-sec 3

