## 04b_contrast_covariance.R -- subtype contrast allowing for the covariance
## contributed by variants shared between the Graves' disease and autoimmune
## hypothyroidism instrument sets.
##
## The primary contrast in 04_sens.R treats the two IVW estimates as
## independent. They share the outcome GWAS and, in part, the FinnGen control
## group, so this script subtracts twice the covariance that can be estimated
## from shared variants. The shared-control component cannot be identified from
## summary data and is therefore not captured; the result is a partial
## sensitivity analysis, not a replacement for the primary test.

source("00_config.R")
suppressMessages({ library(data.table); library(TwoSampleMR) })
ln2 <- log(2)

harm <- function(ex, oc) {
  d <- fread(sprintf("data/mr/h/%s_%s.tsv", ex, oc))
  exd <- d[, c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure",
               "other_allele.exposure", "eaf.exposure", "pval.exposure",
               "exposure", "id.exposure"), with = FALSE]
  ocd <- d[, c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome",
               "other_allele.outcome", "eaf.outcome", "outcome", "id.outcome"),
           with = FALSE]
  h <- as.data.table(suppressMessages(harmonise_data(exd, ocd, action = 2)))
  h[h$mr_keep, ]
}

## IVW with multiplicative random effects; returns the weight sum and the
## inflation factor needed for the covariance
ivw <- function(h) {
  w <- h$beta.exposure^2 / h$se.outcome^2
  r <- h$beta.outcome / h$beta.exposure
  b <- sum(w * r) / sum(w)
  W <- sum(w)
  q <- sum(w * (r - b)^2)
  infl <- max(1, q / (nrow(h) - 1))
  list(b = b, v = infl / W, W = W, infl = infl)
}

ocs <- sub("^GD_", "", sub("\\.tsv$", "",
          basename(list.files("data/mr/h", pattern = "^GD_"))))
res <- list()
for (oc in ocs) {
  if (!file.exists(sprintf("data/mr/h/AIH_%s.tsv", oc))) next
  ha <- harm("GD", oc); hb <- harm("AIH", oc)
  if (nrow(ha) < 3 || nrow(hb) < 3) next
  a <- ivw(ha); b <- ivw(hb)
  m <- merge(ha[, .(SNP, bx = beta.exposure, so = se.outcome)],
             hb[, .(SNP, by = beta.exposure)], by = "SNP")
  cov_fe <- if (nrow(m)) sum(m$bx * m$by / m$so^2) / (a$W * b$W) else 0
  cov <- cov_fe * sqrt(a$infl * b$infl)
  vdiff <- a$v + b$v - 2 * cov
  diff <- a$b - b$b
  res[[oc]] <- data.table(
    outcome = oc,
    diff = diff * ln2,
    se_ind = sqrt(a$v + b$v) * ln2,
    se_cov = sqrt(vdiff) * ln2,
    p_ind = 2 * pnorm(-abs(diff / sqrt(a$v + b$v))),
    p_cov = 2 * pnorm(-abs(diff / sqrt(vdiff))),
    n_shared = nrow(m))
}
R <- rbindlist(res)
R[, q_ind := p.adjust(p_ind, FDR_METHOD)]
R[, q_cov := p.adjust(p_cov, FDR_METHOD)]
setorder(R, p_ind)
print(R[, .(outcome, diff = round(diff, 4), se_ind = round(se_ind, 4),
            se_cov = round(se_cov, 4), p_ind = signif(p_ind, 3),
            p_cov = signif(p_cov, 3), q_ind = signif(q_ind, 3),
            q_cov = signif(q_cov, 3), n_shared)])
fwrite(R, "results/contrast_covariance.tsv", sep = "\t")
