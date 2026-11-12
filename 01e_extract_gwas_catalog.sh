#!/bin/bash
# 01e_extract_gwas_catalog.sh -- extract ferritin and transferrin saturation at
# the hg38 instrument positions. Run after 01a_keys.R.
set -eu
cd "$(dirname "$0")"
test -s data/mr/keys_hg38.txt
python3 python/extract_positions.py data/raw/GCST90270865.tsv.gz     data/mr/keys_hg38.txt ferritin data/mr/out_ferritin.tsv
python3 python/extract_positions.py data/raw/GCST90476445.tsv.gz     data/mr/keys_hg38.txt TSAT data/mr/out_TSAT.tsv
for f in data/mr/out_ferritin.tsv data/mr/out_TSAT.tsv; do
  n=$(($(wc -l < "$f") - 1))
  echo "$f: $n variants"
  [ "$n" -ge 3 ] || { echo "ERROR: too few variants in $f" >&2; exit 1; }
done
