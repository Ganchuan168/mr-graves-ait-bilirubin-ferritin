# Analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24158221.svg)](https://doi.org/10.5281/zenodo.24158221)

Divergent bilirubin and concordant direction of ferritin effect in Graves'
disease and autoimmune hypothyroidism: a subtype-resolved two-sample Mendelian
randomization study

This folder contains the R and Python code used for the analyses reported in the
manuscript. All inputs are publicly available summary statistics; no
individual-level data are used.

---

## 1. Software

| Component | Version used |
|---|---|
| R | 4.6.1 |
| data.table | 1.18.6.1 |
| TwoSampleMR | 0.7.11 |
| ieugwasr | 1.2.0 |
| MVMR | 0.4.8 |
| sandwich | 3.1-3 |
| ggplot2 | 4.0.3 |
| cowplot | 1.2.0 |
| PLINK | 1.9 |
| Python | 3.12.3 (numpy, pandas, matplotlib 3.11.2) |
| curl, bgzip, awk | streaming download and extraction |

`sessionInfo.txt` records the exact R session used for the reported analyses. The
file produced by the final reported run is included at the package root, so the
loaded versions can be checked without rerunning anything. It is
written by `13_session_info.R`, which loads every analysis package in one process
before calling `sessionInfo()`; a bare `Rscript -e` call would record only base R.
`renv.lock` pins R 4.6.1, the seven direct dependencies and their full dependency
closure (91 packages). `package_versions.tsv` lists the same seven direct
dependencies in plain text: data.table, TwoSampleMR, ieugwasr, MVMR, sandwich,
ggplot2 and cowplot. CRAN packages are recorded by version. Three direct
dependencies (TwoSampleMR, ieugwasr, MVMR) and three packages that TwoSampleMR
imports but this pipeline never calls (RadialMR, MR-PRESSO, MRMix) are distributed
through GitHub and are recorded by repository, release tag or commit and commit
SHA. RadialMR, MR-PRESSO and MRMix are present only because TwoSampleMR lists them
under `Imports`; no MR-PRESSO, RadialMR or MRMix analysis was run. `renv::restore()`
recreates the library from the lockfile. The frozen inputs and `14_check_results.R`
verify the numerical results independently of the package manager.
The pipeline was run on Linux (x86_64); it requires network access to FinnGen,
the Pan-UK Biobank S3 bucket and the GWAS Catalog FTP site. A complete run takes
about one to two hours and less than 8 GB of RAM.

## 2. Directory layout expected at run time

```
.
├── 00_config.R … 14_check_results.R      scripts (this folder)
├── python/                                helper scripts
├── frozen/                                harmonised inputs for off-line reproduction
├── data/raw/                              downloaded summary statistics
├── data/mr/                               instruments, extracted and harmonised files
├── results/                               analysis output (TSV)
└── figures/                               generated figures
```

All scripts are run from the package root. `00_config.R` creates the
directories and holds every tunable parameter:

| Object | Value | Meaning |
|---|---|---|
| `LD_REF` | `ref/1kg_eur/EUR` or `$AITD_LD_REF` | PLINK prefix of the 1000 Genomes phase 3 European panel (hg19) |
| `PLINK_BIN` | `plink1.9` or `$AITD_PLINK` | PLINK executable |
| `SEED` | 20261111 | seed for every randomised step |
| `P_PREFILTER` | 1e-5 | streaming pre-filter on the FinnGen p value column |
| `P_INSTRUMENT` | 5e-8 | genome-wide significance threshold |
| `CLUMP_KB`, `CLUMP_R2` | 10000, 0.001 | LD clumping window and threshold |
| `HLA_START`, `HLA_END` | 25e6, 35e6 | extended MHC, chr6, GRCh37 |
| `TRIM_K` | 3 | absolute studentised residual for outlier trimming |
| `FDR_METHOD` | `BH` | Benjamini–Hochberg |
| `NON_OUTCOME_TRAITS` | `CREAT` | serum creatinine is a mediator, excluded from the outcome loop |
| `N_BOOT` | 20000 | parametric bootstrap draws |
| `AITD_CASES`, `AITD_CONTROLS` | see file | released FinnGen R12 counts (Table 1 only) |

Point `AITD_LD_REF` at the `.bed/.bim/.fam` prefix of the 1000 Genomes phase 3
European panel (hg19, about 489 unrelated individuals).

## 3. Input data

### Exposures (FinnGen release R12)

| Exposure | Endpoint | File written by `01d_download_public.sh` |
|---|---|---|
| Graves' disease | `E4_GRAVES_STRICT` | `data/raw/fg_GRAVES.gz` |
| Autoimmune hypothyroidism | `E4_HYTHY_AI_STRICT` | `data/raw/fg_HYTHY_AI.gz` |
| Autoimmune hyperthyroidism | `AUTOIMMUNE_HYPERTHYROIDISM` | `data/raw/fg_AIHYPERTHY.gz` |

### Outcomes

Fourteen blood-cell and biochemistry traits and C-reactive protein come from the
Pan-UK Biobank flat files
(`https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files`), extracted at
the instrument positions in hg19.

| Trait | Pan-UKB field | File |
|---|---|---|
| Reticulocyte count | 30250 | `data/mr/out_retic_count.tsv` |
| Reticulocyte percentage | 30240 | `data/mr/out_retic_pct.tsv` |
| Immature reticulocyte fraction | 30280 | `data/mr/out_IRF.tsv` |
| Mean reticulocyte volume | 30260 | `data/mr/out_MRV.tsv` |
| Red cell distribution width | 30070 | `data/mr/out_RDW.tsv` |
| Mean corpuscular volume | 30040 | `data/mr/out_MCV.tsv` |
| Mean corpuscular haemoglobin | 30050 | `data/mr/out_MCH.tsv` |
| MCH concentration | 30060 | `data/mr/out_MCHC.tsv` |
| Haemoglobin | 30020 | `data/mr/out_Hb.tsv` |
| Red blood cell count | 30010 | `data/mr/out_RBC.tsv` |
| Haematocrit | 30030 | `data/mr/out_Hct.tsv` |
| Total bilirubin | 30840 | `data/mr/out_TBIL.tsv` |
| Direct bilirubin | 30660 | `data/mr/out_DBIL.tsv` |
| C-reactive protein | 30710 | `data/mr/out_CRP.tsv` |
| Serum creatinine | 30700 | `data/mr/out_CREAT.tsv` |

Ferritin and transferrin saturation come from the GWAS Catalog and are
extracted at the hg38 instrument positions.

| Trait | Accession | File |
|---|---|---|
| Ferritin | GCST90270865 | `data/mr/out_ferritin.tsv` |
| Transferrin saturation | GCST90476445 | `data/mr/out_TSAT.tsv` |

Every extracted file has the columns `trait, chr, pos, ea, nea, beta, se, eaf, p`,
with alleles oriented to the effect allele.

## 4. Run order

The whole pipeline runs with:

```bash
export AITD_LD_REF=/path/to/1kg_eur/EUR     # PLINK prefix, hg19
bash run_all.sh
```

The individual steps, in dependency order, are:

```bash
bash 01d_download_public.sh     # raw FinnGen and GWAS Catalog files
Rscript 01_instruments.R        # clumping and instrument selection
Rscript 01a_keys.R              # hg19 and hg38 position lists
bash 01b_fetch_panukb.sh        # Pan-UKB traits at instrument positions
bash 01c_fetch_cis.sh           # the four cis instruments in every trait
bash 01f_cis_aitd.sh            # FinnGen effects at the same four cis variants
bash 01e_extract_gwas_catalog.sh# ferritin and transferrin saturation

Rscript 02_build.R              # merge exposure and outcome files
Rscript 03_mr.R                 # IVW, MR-Egger, weighted median, weighted mode
Rscript 04_sens.R               # FDR, minimum detectable effect, subtype contrast
Rscript 04b_contrast_covariance.R  # contrast with shared-variant covariance
Rscript 05_mvmr.R               # AITD adjusted for CRP
Rscript 06_robust.R             # HLA exclusion, outlier-trimmed IVW
Rscript 09_mvmr3.R              # AITD + CRP + creatinine, three-exposure MVMR
Rscript 10_mediation.R          # two-step MR through CRP and creatinine
Rscript 10b_mediation_bootstrap.R  # bootstrap intervals for proportions

Rscript 07_export.R             # manuscript tables and outcome labels
Rscript 12_variant_table.R      # harmonised variant-level table
Rscript 08_figures.R            # Figures 2 and 3
Rscript 11_figure1_jtm.R        # Figure 1 in the journal layout
python3 python/12_figure4_pathway.py   # Figure 4
Rscript 13_session_info.R       # sessionInfo with all analysis packages loaded
Rscript 14_check_results.R      # verify every reported estimate against expected_*.tsv
```

`07_export.R` must precede `12_variant_table.R`, which reads the outcome-label
file it writes.

### Off-line reproduction from frozen inputs

`frozen/` contains the 48 harmonised exposure–outcome files and the small
mediator and cis extracts. This reproduces every estimate without FinnGen,
Pan-UK Biobank, the GWAS Catalog or an LD reference panel:

```bash
bash run_frozen.sh
```

The script copies `frozen/` into `data/`, runs the analysis, and finishes with
`14_check_results.R`, which compares the 48 primary estimates, the 48 HLA/outlier
sensitivity estimates, 96 heterogeneity tests, 48 Egger intercepts, 20 MVMR
models, 40 mediation rows, 40 bootstrap rows and 16 covariance-adjusted
contrasts with the `expected_*.tsv` files and stops on any mismatch.

`01d` must precede `01_instruments.R` because the FinnGen files are its input;
`01a_keys.R` must precede `01e` because the ferritin and transferrin saturation
extraction uses the hg38 key list.

## 5. What each script does

| Script | Purpose |
|---|---|
| `01_instruments.R` | Reads FinnGen genome-wide summary statistics, keeps variants with P < 1e-5, restricts to variants in the LD panel, clumps at 10 Mb and r² < 0.001, computes F statistics, writes `data/mr/instr_*.tsv` |
| `01a_keys.R` | Writes hg19 and hg38 position lists for the extraction scripts |
| `01b_fetch_panukb.sh`, `stream.sh` | Streams each Pan-UKB flat file through `curl \| bgzip -dc \| awk` and keeps the requested positions; fails on an empty extraction |
| `01c_fetch_cis.sh` | Fetches the three CRP-locus variants (chr1:159.64–159.77 Mb, hg19) and the UMOD variant (chr16:20.39 Mb) in every trait |
| `01f_cis_aitd.sh` | Reads the same four variants by rsID from the FinnGen files and writes `data/mr/cis_aitd.tsv`, the exposure side of the cis block |
| `01d_download_public.sh` | Downloads the FinnGen endpoints and the two GWAS Catalog files, with gzip integrity checks |
| `01e_extract_gwas_catalog.sh` | Extracts ferritin and transferrin saturation at the hg38 instrument positions and checks the row counts |
| `02_build.R` | Maps exposure variants to hg19, merges with each outcome file, writes harmonisation-ready tables to `data/mr/h/` |
| `03_mr.R` | Harmonises each pair (`harmonise_data(action = 2)`) and fits IVW (multiplicative random effects), MR-Egger, weighted median and weighted mode; heterogeneity, pleiotropy, leave-one-out and directionality tests |
| `04_sens.R` | Benjamini–Hochberg adjustment within exposure, minimum detectable effect at 80% power, and the subtype contrast (Graves' disease minus autoimmune hypothyroidism) |
| `04b_contrast_covariance.R` | Repeats the contrast after subtracting twice the covariance contributed by variants shared between the two instrument sets; a partial adjustment that does not capture shared-control correlation |
| `05_mvmr.R` | Two-exposure multivariable MR of AITD liability adjusted for CRP, with `MVMR::strength_mvmr` conditional F statistics |
| `06_robust.R` | Repeats the IVW estimate after excluding HLA-region variants and after trimming variants with studentised residuals above 3 |
| `09_mvmr3.R` | Three-exposure multivariable MR (AITD, CRP, creatinine); combines instrument and cis blocks and reports conditional F statistics |
| `10_mediation.R` | Two-step MR: AITD to mediator, mediator to outcome using cis instruments only, product-of-coefficients indirect effect with delta-method interval |
| `10b_mediation_bootstrap.R` | Parametric bootstrap (20,000 draws) for the mediated proportions |
| `12_variant_table.R` | Exports `results/harmonised_variants.tsv`, the variant-level table that underlies Additional file 4; HLA flags use hg19 coordinates and each variant carries an approximate R² |
| `13_session_info.R` | Loads all analysis packages and writes a genuine `sessionInfo.txt` |
| `14_check_results.R` | Compares every primary estimate, sensitivity estimate, heterogeneity test, Egger intercept, MVMR model and mediation row with the `expected_*.tsv` files; fails loudly on any difference |
| `07_export.R` | Converts log-odds to per-doubling effects (multiplication by ln 2) and writes the manuscript tables |
| `08_figures.R` | Forest plot of subtype contrasts and variant-level scatter plots (Figures 2 and 3) |
| `11_figure1_jtm.R` | Figure 1: forest plot of all estimates by outcome domain |
| `12_figure4_pathway.py` | Figure 4: pathway diagram and mediated proportions |

## 6. Output to manuscript mapping

| Manuscript item | File in `results/` |
|---|---|
| Table 1 | `table1_instruments.tsv` |
| Table 2 | `table2_primary.tsv` |
| Table 3 (subtype contrasts) | `table3_contrast.tsv` |
| Table 4 (sensitivity analyses) | `table4_sensitivity.tsv`, `robustness.tsv` |
| Table 5 (multivariable and two-step MR) | `mvmr3.tsv`, `mediation.tsv` |
| Additional file 1, Tables S1–S5 | `primary_summary.tsv`, `mvmr.tsv`, `leaveoneout.tsv`, `mediation.tsv`, `robustness.tsv` |
| Additional file 1, Table S6 | `contrast_covariance.tsv` |
| Additional file 1, Table S7 | `mediation_bootstrap.tsv` |
| Additional file 1, Table S8 | built from the field lists in section 3 |
| Additional file 4 | `harmonised_variants.tsv`, written by `12_variant_table.R` |
| Figure 1 | `figures/fig1_forest_jtm.*` |
| Figure 2 | `figures/fig2_contrast.*` |
| Figure 3 | `figures/fig3_scatter.*` |
| Figure 4 | `figures/fig4_pathway.png`, `figures/fig4_pathway.pdf` |

The file names now match the manuscript table numbers.

## 7. Statistical notes and known behaviour

- All effects are expressed per doubling of genetic liability (the IVW
  coefficient multiplied by ln 2). Outcomes are inverse-rank-normalised, so
  coefficients are in standard deviation units.
- `data/mr/h/*.tsv` are merged but **not** harmonised: `03_mr.R` and the
  sensitivity scripts call `harmonise_data(action = 2)` at read time. This
  removes variants with allele mismatch or missing outcome data and palindromic
  variants whose allele frequency lies between 0.42 and 0.58; palindromic
  variants at extreme frequencies are retained and aligned. Additional file 4
  records the keep indicator and the reason for every excluded variant.
- The subtype contrast assumes the two IVW estimates are independent. They are
  not: they share the Pan-UK Biobank outcome data and part of the FinnGen
  control group. `04b_contrast_covariance.R` subtracts the covariance
  attributable to shared variants (one variant, rs11571297) but cannot estimate
  the shared-control component; the primary test therefore remains the
  independence test, described as conservative under positive correlation.
- The two-step mediation uses disjoint variant sets in the two steps, which
  prevents variant overlap, but the mediator and outcome GWAS both come from the
  Pan-UK Biobank and therefore overlap in participants. The second step is
  biased towards the observational mediator–outcome association, and the
  mediation results are reported as exploratory.
- Mediated proportions are ratios with an imprecisely estimated denominator.
  The delta-method interval can extend outside [0, 1]; the bootstrap interval in
  `10b_mediation_bootstrap.R` makes this explicit, and the proportions are
  interpreted descriptively.
- Randomised steps use `SEED` from `00_config.R`. `13_session_info.R` writes
  `sessionInfo.txt` after loading the analysis packages.
- The scripts were written for a Linux environment; paths use `/` and the
  streaming step needs `curl`, `bgzip` and `awk`.

## 8. Data and code availability

The code is released under the MIT licence (see `LICENSE`). The complete code,
the frozen inputs and the harmonised variant-level tables are supplied as
Additional files 3 and 4 and form part of the article. The complete reproducible
R analysis code package has been assigned a reserved Zenodo DOI
(10.5281/zenodo.24158221), which will resolve when the record is made public upon
publication. No request to the authors is needed. `CITATION.cff` and `.zenodo.json`
carry the title, author order, ORCID iDs and licence for the Zenodo record. Source GWAS data remain subject to the terms of FinnGen,
Pan-UK Biobank and the GWAS Catalog.

## 9. Contact

Meng Xia, School of Basic Medicine, Guangxi University of Chinese Medicine,
Nanning, China. E-mail: xiam@gxtcmu.edu.cn
