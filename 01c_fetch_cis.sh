#!/bin/bash
# 01c_fetch_cis.sh -- four cis instruments (three at the CRP locus, one at UMOD)
# in every outcome, used by the multivariable and two-step analyses.
set -eu
cd "$(dirname "$0")"
cat > data/mr/keys_cis.txt <<'KEYS'
1:159642349
1:159723120
1:159773350
16:20392332
KEYS

B="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"
for p in "30250 retic_count cont" "30240 retic_pct cont" "30020 Hb cont" "30030 Hct cont" \
         "30010 RBC cont" "30040 MCV cont" "30260 MRV cont" "30070 RDW cont" \
         "30840 TBIL bio" "30660 DBIL bio" "30710 CRP bio" "30700 CREAT bio"; do
  set -- $p
  python3 python/tabix_fetch.py "$B/continuous-$1-both_sexes-irnt.tsv.bgz" \
      data/mr/keys_cis.txt "data/mr/cis_$2.tsv" "$2" \
    || python3 python/tabix_fetch.py "$B/biomarkers-$1-both_sexes-irnt.tsv.bgz" \
      data/mr/keys_cis.txt "data/mr/cis_$2.tsv" "$2"
done
echo "cis extraction finished"
