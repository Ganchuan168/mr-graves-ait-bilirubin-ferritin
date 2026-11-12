## 04_sens.R -- sensitivity, subtype contrast, FDR, HLA exclusion, MDE
source("00_config.R")
suppressMessages({ library(data.table); library(TwoSampleMR) })

lab <- c(retic_count="Reticulocyte count", retic_pct="Reticulocyte percentage",
  IRF="Immature reticulocyte fraction", MRV="Mean reticulocyte volume",
  RDW="Red cell distribution width", MCV="Mean corpuscular volume",
  MCH="Mean corpuscular haemoglobin", MCHC="MCHC", Hb="Haemoglobin",
  RBC="Red blood cell count", Hct="Haematocrit", TBIL="Total bilirubin",
  DBIL="Direct bilirubin", ferritin="Ferritin", TSAT="Transferrin saturation",
  CRP="C-reactive protein")
exlab <- c(GD="Graves disease", AIH="Autoimmune hypothyroidism", AIT="Autoimmune hyperthyroidism")

files <- list.files("data/mr/h", pattern = "\\.tsv$", full.names = TRUE)
out <- list()
for (f in files) {
  d <- fread(f); if (nrow(d) < 5) next
  ex <- d$id.exposure[1]; oc <- d$id.outcome[1]
  exd <- d[, c("SNP","beta.exposure","se.exposure","effect_allele.exposure",
               "other_allele.exposure","eaf.exposure","pval.exposure","exposure","id.exposure"), with=FALSE]
  ocd <- d[, c("SNP","beta.outcome","se.outcome","effect_allele.outcome",
               "other_allele.outcome","eaf.outcome","outcome","id.outcome"), with=FALSE]
  h <- harmonise_data(exd, ocd, action = 2); h <- h[h$mr_keep, ]
  if (nrow(h) < 5) next
  ## mr_ivw always fits the multiplicative random-effects model in
  ## TwoSampleMR >= 0.5; no method string is needed (or read) here.
  ivw_all <- mr_ivw(h$beta.exposure, h$beta.outcome, h$se.exposure, h$se.outcome)
  fstat <- (h$beta.exposure / h$se.exposure)^2
  out[[length(out)+1]] <- data.table(
    exposure = ex, outcome = oc, n = nrow(h),
    b = ivw_all$b, se = ivw_all$se, p = ivw_all$pval,
    mde80 = (qnorm(0.975) + qnorm(0.8)) * ivw_all$se,
    meanF = mean(fstat), minF = min(fstat))
}
r <- rbindlist(out)
r[, q := p.adjust(p, FDR_METHOD), by = exposure]
fwrite(r, "results/primary_summary.tsv", sep = "\t")
print(r[order(exposure, outcome), .(exposure, outcome, n, b = round(b,4), se = round(se,4),
                                    p = signif(p,3), q = signif(q,3), mde80 = round(mde80,3))], nrows = 100)

## subtype contrast for shared outcomes
w <- dcast(r, outcome ~ exposure, value.var = c("b","se"))
w <- w[!is.na(b_GD) & !is.na(b_AIH)]
w[, diff := b_GD - b_AIH]
w[, se_diff := sqrt(se_GD^2 + se_AIH^2)]
w[, z := diff / se_diff]
w[, p_diff := 2 * pnorm(-abs(z))]
fwrite(w, "results/subtype_contrast.tsv", sep = "\t")
cat("\nSubtype contrast (GD vs AIH):\n")
print(w[, .(outcome, b_GD = round(b_GD,4), b_AIH = round(b_AIH,4),
            diff = round(diff,4), p_diff = signif(p_diff,3))])
