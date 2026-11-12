## 07_export.R -- manuscript tables (per-doubling scaling, ln2)
source("00_config.R")
suppressMessages(library(data.table))
ln2 <- log(2)

lab <- data.table(
  outcome = c("retic_count","retic_pct","IRF","Hb","Hct","RBC","MCV","MCH","MCHC","RDW","MRV",
         "TBIL","DBIL","ferritin","TSAT","CRP"),
  name = c("Reticulocyte count","Reticulocyte percentage","Immature reticulocyte fraction",
           "Haemoglobin","Haematocrit","Red blood cell count","Mean corpuscular volume",
           "Mean corpuscular haemoglobin","Mean corpuscular haemoglobin concentration",
           "Red cell distribution width","Mean reticulocyte volume",
           "Total bilirubin","Direct bilirubin","Ferritin","Transferrin saturation",
           "C-reactive protein"),
  domain = c(rep("Erythropoietic output",3), rep("Red cell indices",8),
             rep("Haemolysis / heme catabolism",2), rep("Iron status",2), "Inflammation"))
lab[, dom := factor(domain, levels = c("Erythropoietic output","Red cell indices",
                                       "Haemolysis / heme catabolism","Iron status","Inflammation"))]
setorder(lab, dom)
fwrite(lab[, .(outcome, name, domain, order = as.integer(dom))], "results/outcome_labels.tsv", sep = "\t")

## ---- Table 1: instruments ----
t1 <- rbindlist(lapply(c("GD","AIH","AIT"), function(x) {
  d <- fread(sprintf("data/mr/instr_%s.tsv", x))
  data.table(exposure = x, n_snp = nrow(d), mean_F = mean(d$fstat),
             median_F = median(d$fstat), min_F = min(d$fstat),
             n_F20 = sum(d$fstat > 20), n_F30 = sum(d$fstat > 30),
             n_F50 = sum(d$fstat > 50))
}))
t1[, exposure := c(GD = "Graves disease", AIH = "Autoimmune hypothyroidism",
                   AIT = "Autoimmune hyperthyroidism")[exposure]]
t1[, cases := AITD_CASES[exposure]]
t1[, controls := AITD_CONTROLS[exposure]]
fwrite(t1, "results/table1_instruments.tsv", sep = "\t")

## ---- Table 2: primary IVW per doubling ----
p <- fread("results/primary_summary.tsv")
p <- merge(p, lab[, .(outcome, name, domain, dom)], by = "outcome", all.x = TRUE)
p[, `:=`(b_d = b * ln2, lo = (b - 1.96 * se) * ln2, hi = (b + 1.96 * se) * ln2,
         mde_d = mde80 * ln2)]
p[, exposure := c(GD = "Graves disease", AIH = "Autoimmune hypothyroidism",
                  AIT = "Autoimmune hyperthyroidism")[exposure]]
setorder(p, dom, outcome, exposure)
fwrite(p[, .(domain, outcome = name, exposure, n, b_d, lo, hi, p, q, mde_d)],
       "results/table2_primary.tsv", sep = "\t")

## ---- Table 3: sensitivity ----
m <- fread("results/mr_estimates.tsv")
map <- c("Inverse variance weighted (multiplicative random effects)" = "IVW",
         "MR Egger" = "Egger", "Weighted median" = "WMed", "Weighted mode" = "WMode")
m <- m[method %in% names(map)]
m[, mm := map[method]]
b <- dcast(m, id.exposure + id.outcome ~ mm, value.var = "b")
pv <- dcast(m, id.exposure + id.outcome ~ mm, value.var = "pval")
setnames(b, 3:6, paste0("b", names(b)[3:6]))
setnames(pv, 3:6, paste0("p", names(pv)[3:6]))
setnames(b, c("id.exposure","id.outcome"), c("exposure","outcome"))
setnames(pv, c("id.exposure","id.outcome"), c("exposure","outcome"))
h <- fread("results/heterogeneity.tsv")[
  method == "Inverse variance weighted (multiplicative random effects)",
  .(exposure = id.exposure, outcome = id.outcome, Q = Q, Q_df = Q_df, Q_pval)]
pl <- fread("results/pleiotropy.tsv")[
  , .(exposure = id.exposure, outcome = id.outcome, int = egger_intercept, int_p = pval)]
s <- Reduce(function(a, c) merge(a, c, by = c("exposure","outcome"), all = TRUE), list(b, pv, h, pl))

rb <- fread("results/robustness.tsv")
s <- merge(s, rb[, .(exposure, outcome, b_nohla, se_nohla, p_nohla, b_trim, p_trim, n_trimmed)],
           by = c("exposure","outcome"), all.x = TRUE)

s[, `:=`(bIVW = bIVW * ln2, bEgger = bEgger * ln2, bWMed = bWMed * ln2,
         bWMode = bWMode * ln2, b_nohla = b_nohla * ln2, b_trim = b_trim * ln2)]
s <- merge(s, lab[, .(outcome, name, dom)], by = "outcome", all.x = TRUE)
s[, exposure := c(GD = "Graves disease", AIH = "Autoimmune hypothyroidism",
                  AIT = "Autoimmune hyperthyroidism")[exposure]]
setorder(s, dom, outcome, exposure)
fwrite(s, "results/table4_sensitivity.tsv", sep = "\t")

## ---- Table 3: subtype contrast (all and HLA-excluded) ----
w <- fread("results/subtype_contrast.tsv")
w <- merge(w, lab, by = "outcome", all.x = TRUE)
w[, `:=`(gd = b_GD * ln2, gd_lo = (b_GD - 1.96 * se_GD) * ln2, gd_hi = (b_GD + 1.96 * se_GD) * ln2,
         ai = b_AIH * ln2, ai_lo = (b_AIH - 1.96 * se_AIH) * ln2, ai_hi = (b_AIH + 1.96 * se_AIH) * ln2,
         diff = diff * ln2,
         diff_lo = (b_GD - b_AIH - 1.96 * se_diff) * ln2,
         diff_hi = (b_GD - b_AIH + 1.96 * se_diff) * ln2)]

rbg <- rb[exposure == "GD", .(outcome, gd_nh = b_nohla * ln2, seg = se_nohla * ln2)]
rba <- rb[exposure == "AIH", .(outcome, ai_nh = b_nohla * ln2, sea = se_nohla * ln2)]
w <- merge(w, rbg, by = "outcome", all.x = TRUE)
w <- merge(w, rba, by = "outcome", all.x = TRUE)
w[, `:=`(diff_nohla = gd_nh - ai_nh, se_nh = sqrt(seg^2 + sea^2))]
w[, p_nohla := 2 * pnorm(-abs(diff_nohla / se_nh))]
w[, `:=`(diff_nohla_lo = diff_nohla - 1.96 * se_nh, diff_nohla_hi = diff_nohla + 1.96 * se_nh)]
w[, q_diff := p.adjust(p_diff, method = FDR_METHOD)]
setorder(w, dom, outcome)
fwrite(w[, .(domain, outcome = name, gd, gd_lo, gd_hi, ai, ai_lo, ai_hi,
             diff, diff_lo, diff_hi, p_diff, q_diff, diff_nohla, diff_nohla_lo, diff_nohla_hi, p_nohla)],
       "results/table3_contrast.tsv", sep = "\t")

cat("exported\n")
print(t1)
print(w[, .(outcome = name, gd = round(gd,4), ai = round(ai,4), diff = round(diff,4),
            p_diff = signif(p_diff,3), p_nohla = signif(p_nohla,3))])
