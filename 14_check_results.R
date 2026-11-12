## 14_check_results.R -- verify that a re-run reproduces the primary estimates,
## the sensitivity analyses, the multivariable models and the mediation outputs.
source("00_config.R")
suppressMessages(library(data.table))

cmp <- function(exp_file, got_file, keys, cols, label) {
  e <- fread(exp_file); g <- fread(got_file)
  m <- merge(e, g, by = keys, suffixes = c(".exp", ".got"))
  if (nrow(m) != nrow(e)) stop(label, ": row count differs after merge")
  worst <- 0
  for (cc in cols) {
    a <- m[[paste0(cc, ".exp")]]; b <- m[[paste0(cc, ".got")]]
    ok <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) &
          abs(a - b) <= 1e-6 + 1e-6 * abs(a))
    if (!all(ok)) stop(label, ": ", cc, " does not reproduce")
    d <- suppressWarnings(max(abs(a - b), na.rm = TRUE))
    if (is.finite(d)) worst <- max(worst, d)
  }
  cat(sprintf("%-22s %3d rows; max |difference| = %.3g\n", label, nrow(m), worst))
}

cmp("expected_results.tsv", "results/primary_summary.tsv",
    c("exposure", "outcome"), c("b", "se", "p", "n"), "primary IVW")
cmp("expected_robustness.tsv", "results/robustness.tsv",
    c("exposure", "outcome"),
    c("n", "n_hla", "b", "p", "b_nohla", "p_nohla", "b_trim", "p_trim", "n_trimmed"),
    "HLA / outlier")
cmp("expected_heterogeneity.tsv", "results/heterogeneity.tsv",
    c("exposure", "outcome", "method"), c("Q", "Q_df", "Q_pval"), "heterogeneity")
cmp("expected_pleiotropy.tsv", "results/pleiotropy.tsv",
    c("exposure", "outcome"), c("egger_intercept", "se", "pval"), "Egger intercept")
cmp("expected_mvmr3.tsv", "results/mvmr3.tsv",
    c("exposure", "outcome"),
    c("n", "n_cis", "bAITD", "seAITD", "pAITD", "bCRP", "pCRP", "bKID", "pKID",
      "condF_AITD", "condF_CRP", "condF_KID"), "three-exposure MVMR")
cmp("expected_mediation.tsv", "results/mediation.tsv",
    c("exposure", "outcome", "mediator"),
    c("b1", "se1", "b2", "se2", "ind", "se_ind", "prop"), "two-step mediation")
cmp("expected_bootstrap.tsv", "results/mediation_bootstrap.tsv",
    c("exposure", "outcome", "mediator"),
    c("prop_boot_lo", "prop_boot_hi", "prop_boot_median"), "mediation bootstrap")
cmp("expected_contrast_cov.tsv", "results/contrast_covariance.tsv",
    c("outcome"), c("diff", "se_cov", "p_cov", "q_cov"), "contrast covariance")

cat("All reported primary and sensitivity estimates reproduced.\n")
