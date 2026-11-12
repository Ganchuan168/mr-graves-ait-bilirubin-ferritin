#!/bin/bash
# run_all.sh -- execute the pipeline in dependency order.
# The LD reference must be set first:
#   export AITD_LD_REF=/path/to/1kg_eur/EUR
set -eu
cd "$(dirname "$0")"

bash 01d_download_public.sh          # raw FinnGen and GWAS Catalog files
Rscript 01_instruments.R             # clumping and instrument selection
Rscript 01a_keys.R                   # hg19/hg38 position lists
bash 01b_fetch_panukb.sh             # Pan-UKB outcomes at instrument positions
bash 01c_fetch_cis.sh                # cis instruments in every trait
bash 01f_cis_aitd.sh                 # FinnGen effects at the four cis variants
bash 01e_extract_gwas_catalog.sh     # ferritin and transferrin saturation

Rscript 02_build.R                   # merge exposure and outcome files
Rscript 03_mr.R                      # IVW, MR-Egger, weighted median/mode
Rscript 04_sens.R                    # FDR, MDE, subtype contrast
Rscript 04b_contrast_covariance.R    # contrast with shared-variant covariance
Rscript 05_mvmr.R                    # two-exposure MVMR adjusted for CRP
Rscript 06_robust.R                  # HLA exclusion and outlier trimming
Rscript 09_mvmr3.R                   # three-exposure MVMR and two-step mediation
Rscript 10_mediation.R               # product-of-coefficients mediation
Rscript 10b_mediation_bootstrap.R    # bootstrap intervals for proportions

Rscript 07_export.R                  # manuscript tables and outcome labels
Rscript 12_variant_table.R           # harmonised variant-level table
Rscript 08_figures.R                 # Figures 2 and 3
Rscript 11_figure1_jtm.R             # Figure 1
python3 python/12_figure4_pathway.py # Figure 4
Rscript 13_session_info.R            # genuine sessionInfo with all packages loaded
Rscript 14_check_results.R           # compare key estimates with expected_results.tsv
echo "pipeline finished"
