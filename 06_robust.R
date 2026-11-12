## 06_robust.R -- HLA-region exclusion and outlier-trimmed IVW
source("00_config.R")
suppressMessages({ library(data.table); library(TwoSampleMR) })

## hg19 positions are carried in the harmonised files by 02_build.R, so the
## HLA flag does not depend on an external LD panel and cannot mix genome
## builds. The fallback keeps older cached files usable.
bim <- NULL
if (!all(c("chr19", "pos19") %in% names(fread(list.files("data/mr/h",
    pattern = "\\.tsv$", full.names = TRUE)[1], nrows = 1)))) {
  bim <- fread(paste0(LD_REF, ".bim"), header = FALSE)
  setnames(bim, c("chr19", "SNP", "cm", "pos19", "a1", "a2"))
}

ivw <- function(d) {
  if (nrow(d) < 5) return(list(b = NA, se = NA, p = NA, n = nrow(d)))
  r <- mr_ivw(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome)
  list(b = r$b, se = r$se, p = r$p, n = nrow(d))
}

## outlier-trimmed IVW: drop SNPs with |studentized residual| > 3, refit
trim_ivw <- function(d, k = TRIM_K) {
  if (nrow(d) < 6) return(list(b = NA, se = NA, p = NA, n_rm = NA))
  fit <- lm(beta.outcome ~ beta.exposure - 1, data = d, weights = 1 / d$se.outcome^2)
  rs <- rstandard(fit)
  keep <- abs(rs) <= k
  if (sum(keep) == nrow(d)) return(list(b = NA, se = NA, p = NA, n_rm = 0L))
  r <- ivw(d[keep])
  list(b = r$b, se = r$se, p = r$p, n_rm = sum(!keep))
}

files <- list.files("data/mr/h", pattern = "\\.tsv$", full.names = TRUE)
res <- list()
for (f in files) {
  d <- fread(f); if (nrow(d) < 5) next
  ex <- d$id.exposure[1]; oc <- d$id.outcome[1]
  exd <- d[, c("SNP","beta.exposure","se.exposure","effect_allele.exposure",
               "other_allele.exposure","eaf.exposure","pval.exposure","exposure","id.exposure"), with=FALSE]
  ocd <- d[, c("SNP","beta.outcome","se.outcome","effect_allele.outcome",
               "other_allele.outcome","eaf.outcome","outcome","id.outcome"), with=FALSE]
  h <- harmonise_data(exd, ocd, action = 2)
  h <- as.data.table(h[h$mr_keep, ])
  if (nrow(h) < 5) next
  if (!is.null(bim)) {
    h <- merge(h, bim[, .(SNP, chr19, pos19)], by = "SNP", all.x = TRUE)
  } else {
    h <- merge(h, unique(d[, .(SNP, chr19, pos19)]), by = "SNP", all.x = TRUE)
  }
  h[, hla := !is.na(chr19) & chr19 == 6 & pos19 >= HLA_START & pos19 <= HLA_END]

  a <- ivw(h)
  b <- ivw(h[hla == FALSE])
  t <- trim_ivw(h)
  res[[length(res) + 1]] <- data.table(
    exposure = ex, outcome = oc, n = a$n, n_hla = sum(h$hla),
    b = a$b, se = a$se, p = a$p,
    b_nohla = b$b, se_nohla = b$se, p_nohla = b$p, n_nohla = b$n,
    b_trim = t$b, p_trim = t$p, n_trimmed = t$n_rm)
}
R <- rbindlist(res)
fwrite(R, "results/robustness.tsv", sep = "\t")
cat("done", nrow(R), "sets\n")
print(R[order(exposure, outcome), .(exposure, outcome, n, n_hla,
      b = round(b, 4), p = signif(p, 3),
      b_nohla = round(b_nohla, 4), p_nohla = signif(p_nohla, 3),
      b_trim = round(b_trim, 4), p_trim = signif(p_trim, 3))], nrows = 60)
