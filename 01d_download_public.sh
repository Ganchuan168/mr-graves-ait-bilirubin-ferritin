#!/bin/bash
# 01d_download_public.sh -- download FinnGen R12 exposures and the two GWAS
# Catalog traits. Extraction of the GWAS Catalog traits is done by
# 01e_extract_gwas_catalog.sh, which needs the hg38 key list from 01a_keys.R.
set -eu
cd "$(dirname "$0")"
mkdir -p data/raw

dl () {  # $1 url, $2 output
  [ -s "$2" ] && { echo "SKIP $2"; return; }
  curl -fL --retry 5 --retry-delay 3 -o "$2.part" "$1"
  gzip -t "$2.part"
  mv "$2.part" "$2"
  echo "DONE $2 ($(gzip -dc "$2" | wc -l) lines)"
}

F="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
dl "$F/finngen_R12_E4_GRAVES_STRICT.gz"             data/raw/fg_GRAVES.gz
dl "$F/finngen_R12_E4_HYTHY_AI_STRICT.gz"           data/raw/fg_HYTHY_AI.gz
dl "$F/finngen_R12_AUTOIMMUNE_HYPERTHYROIDISM.gz"   data/raw/fg_AIHYPERTHY.gz

E="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
dl "$E/GCST90270001-GCST90271000/GCST90270865/GCST90270865.tsv.gz" data/raw/GCST90270865.tsv.gz
dl "$E/GCST90476001-GCST90477000/GCST90476445/GCST90476445.tsv.gz" data/raw/GCST90476445.tsv.gz
echo "download finished"
