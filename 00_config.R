## 00_config.R -- shared settings for the AITD-anaemia MR pipeline
##
## Run every script from the package root so that data/, results/ and figures/
## are resolved relative to it. The only external path is the LD reference
## panel, which is set here or through the environment variable AITD_LD_REF.

## PLINK binary prefix of the 1000 Genomes phase 3 European panel (hg19),
## used for clumping and for mapping variants to hg19 coordinates.
LD_REF <- Sys.getenv("AITD_LD_REF", "ref/1kg_eur/EUR")

## PLINK executable. Override with AITD_PLINK if plink1.9/plink is elsewhere.
PLINK_BIN <- Sys.getenv("AITD_PLINK",
  unname(Sys.which(c("plink1.9", "plink"))[1]))

## Random seed for every randomised step.
SEED <- 20261111

## Instrument selection
P_PREFILTER   <- 1e-5    # streaming pre-filter on the FinnGen p value column
P_INSTRUMENT  <- 5e-8    # genome-wide significance threshold
CLUMP_KB      <- 10000   # 10 Mb
CLUMP_R2      <- 0.001

## Sensitivity-analysis settings
HLA_START <- 25e6        # chr6:25-35 Mb, GRCh37
HLA_END   <- 35e6
TRIM_K    <- 3           # absolute studentised residual threshold
FDR_METHOD <- "BH"

## Creatinine is a mediator only and is excluded from the outcome loop in
## 02_build.R. CRP is both an outcome and a mediator, so it is not listed here.
NON_OUTCOME_TRAITS <- c("CREAT")

## Parametric bootstrap for mediated proportions
N_BOOT <- 20000

## Released FinnGen R12 case and control counts, used in Table 1 only
AITD_CASES <- c("Graves disease" = 3962,
                "Autoimmune hypothyroidism" = 54752,
                "Autoimmune hyperthyroidism" = 2469)
AITD_CONTROLS <- c("Graves disease" = 496386,
                   "Autoimmune hypothyroidism" = 362639,
                   "Autoimmune hyperthyroidism" = 370637)

for (d in c("data/raw", "data/mr", "data/mr/h", "results", "figures"))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

if (Sys.getenv("AITD_FROZEN") == "" && !file.exists(paste0(LD_REF, ".bim")))
  warning("LD reference not found at ", LD_REF,
          " -- set AITD_LD_REF to the PLINK prefix of the 1000G EUR panel")
