#!/bin/bash
# run_frozen.sh -- reproduce every estimate without downloading the source
# GWAS files. The frozen/ directory contains the harmonised exposure-outcome
# files (frozen/h) and the small mediator/cis extracts (frozen/mr) that the
# pipeline would otherwise rebuild from FinnGen and Pan-UK Biobank.
#
#   bash run_frozen.sh
#
# No LD reference panel is needed in this mode: 02_build.R has already stored
# hg19 coordinates in the harmonised files.
set -eu
cd "$(dirname "$0")"
export AITD_FROZEN=1

[ -d frozen/h ] || { echo "frozen/h is missing" >&2; exit 1; }
mkdir -p data/mr/h results figures
cp frozen/h/*.tsv data/mr/h/
cp frozen/mr/*.tsv data/mr/

Rscript 03_mr.R
Rscript 04_sens.R
Rscript 04b_contrast_covariance.R
Rscript 05_mvmr.R
Rscript 06_robust.R
Rscript 09_mvmr3.R
Rscript 10_mediation.R
Rscript 10b_mediation_bootstrap.R
Rscript 07_export.R
Rscript 12_variant_table.R
Rscript 08_figures.R
Rscript 11_figure1_jtm.R
[ -f python/12_figure4_pathway.py ] && python3 python/12_figure4_pathway.py || true
Rscript 13_session_info.R
Rscript 14_check_results.R
echo "frozen pipeline finished"
