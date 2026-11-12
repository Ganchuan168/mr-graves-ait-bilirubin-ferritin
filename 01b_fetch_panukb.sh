#!/bin/bash
# 01b_fetch_panukb.sh -- extract Pan-UK Biobank summary statistics at the
# instrument and cis-instrument positions (hg19).
set -eu
cd "$(dirname "$0")"
mkdir -p data/mr

while read -r id name layout; do
  [ -z "$id" ] && continue
  case "$id" in \#*) continue ;; esac
  bash stream.sh "$id" "$name" data/mr/keys_hg19.txt "$layout"
done <<'LIST'
30250 retic_count cont
30240 retic_pct   cont
30280 IRF         cont
30260 MRV         cont
30070 RDW         cont
30040 MCV         cont
30050 MCH         cont
30060 MCHC        cont
30020 Hb          cont
30010 RBC         cont
30030 Hct         cont
30840 TBIL        bio
30660 DBIL        bio
30710 CRP         bio
30700 CREAT       bio
LIST
echo "Pan-UK Biobank extraction finished"
