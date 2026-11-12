#!/bin/bash
# stream.sh <field id> <short name> <key file> [cont|bio]
# Streams one Pan-UK Biobank flat file and keeps the requested positions.
set -euo pipefail
field=$1; name=$2; keyfile=$3; layout=${4:-cont}
B="https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files"
if [ "$layout" = "bio" ]; then
  url="$B/biomarkers-${field}-both_sexes-irnt.tsv.bgz"
  AF=14; BE=20; SE=26; PV=32
else
  url="$B/continuous-${field}-both_sexes-irnt.tsv.bgz"
  AF=19; BE=25; SE=31; PV=37
fi
out="data/mr/out_${name}.tsv"
tmp="${out}.tmp"
printf "trait\tchr\tpos\tea\tnea\tbeta\tse\teaf\tp\n" > "$tmp"
curl -fsS --retry 5 --retry-delay 3 "$url" \
 | bgzip -dc \
 | awk -F'\t' -v OFS='\t' -v n="$name" -v af="$AF" -v be="$BE" -v se="$SE" -v pv="$PV" \
       'NR==FNR{k[$1]=1; next}
        FNR==1{next}
        (($1":"$2) in k){ print n,$1,$2,$4,$3,$be,$se,$af,$pv }' \
        "$keyfile" - >> "$tmp"
n=$(($(wc -l < "$tmp") - 1))
[ "$n" -gt 0 ] || { echo "ERROR: no variants extracted for $name" >&2; exit 1; }
mv "$tmp" "$out"
echo "DONE $name $n rows"
