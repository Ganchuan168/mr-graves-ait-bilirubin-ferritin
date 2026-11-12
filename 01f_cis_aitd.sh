#!/bin/bash
# 01f_cis_aitd.sh -- FinnGen R12 effects for the four cis variants used by the
# multivariable and two-step models. 01c_fetch_cis.sh obtains the same four
# variants from Pan-UK Biobank (hg19); this script obtains the exposure side
# from the FinnGen raw files (hg38) and writes them against the Pan-UKB hg19
# key, so data/mr/cis_aitd.tsv is generated rather than assumed.
set -eu
cd "$(dirname "$0")"

## The cis panel is fixed by the design (three CRP-locus variants and one UMOD
## variant). Columns: rsID and Pan-UK Biobank hg19 key.
cat > data/mr/cis_panel.tsv <<'PANEL'
rs1811472	1:159642349
rs6677719	1:159723120
rs61821624	1:159773350
rs77924615	16:20392332
PANEL

printf 'exposure\tSNP\tkey\tea\tnea\tbeta\tse\n' > data/mr/cis_aitd.tsv
cut -f1 data/mr/cis_panel.tsv > data/mr/cis_panel.rsid

extract () {
  local ex=$1 file=$2
  [ -f "$file" ] || { echo "missing $file" >&2; exit 1; }
  ## grep -m4 stops once the four variants have been seen, which avoids
  ## streaming 21 million rows per file through awk.
  zcat "$file" \
    | grep -m4 -w -F -f data/mr/cis_panel.rsid \
    | awk -F'\t' -v OFS='\t' -v ex="$ex" '
        NR==FNR { key[$1]=$2; next }
        {
          split($5, rs, ",");
          print ex, rs[1], key[rs[1]], $4, $3, $9, $10
        }' data/mr/cis_panel.tsv - >> data/mr/cis_aitd.tsv
}

extract GD  data/raw/fg_GRAVES.gz
extract AIH data/raw/fg_HYTHY_AI.gz
[ -f data/raw/fg_AIHYPERTHY.gz ] && extract AIT data/raw/fg_AIHYPERTHY.gz || true

n=$(($(wc -l < data/mr/cis_aitd.tsv) - 1))
[ "$n" -ge 8 ] || { echo "ERROR: only $n cis rows extracted" >&2; exit 1; }
echo "cis_aitd.tsv: $n rows"
