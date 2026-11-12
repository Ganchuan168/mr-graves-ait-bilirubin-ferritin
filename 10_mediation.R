## 10_mediation.R -- two-step MR: AITD -> (CRP | creatinine) -> erythroid outcome
source("00_config.R")
suppressMessages({ library(data.table) })
ln2 <- log(2)
read9 <- function(f) {
  d <- fread(f)
  setnames(d, c("trait", "chr", "pos", "ea", "nea", "beta", "se", "eaf", "p"))
  d[, key := paste(chr, pos, sep = ":")]
  d
}
al <- function(a1, a2, ea, nea, b) {
  out <- rep(NA_real_, length(b))
  i <- ea == a1 & nea == a2; j <- ea == a2 & nea == a1
  out[i] <- b[i]; out[j] <- -b[j]
  out
}

## ---- step 1: AITD -> mediator, AITD instruments only, IVW -----------------
step1 <- function(ex, med) {
  A <- fread(sprintf("data/mr/instr_%s.tsv", ex))[
    , .(SNP, chr19, pos19, ea, nea, beta, se)]
  A[, `:=`(key = paste(chr19, pos19, sep = ":"), a1 = ea, a2 = nea)]
  M <- read9(sprintf("data/mr/out_%s.tsv", ifelse(med == "c", "CRP", "CREAT")))
  z <- merge(A, M[, .(key, ea.m = ea, nea.m = nea, bm = beta, sem = se)], by = "key")
  z[, b_m := al(a1, a2, ea.m, nea.m, bm)]
  z <- z[!is.na(b_m) & !is.na(beta)]
  w <- z$beta^2 / z$sem^2
  b <- sum(w * (z$b_m / z$beta)) / sum(w)
  se <- sqrt(1 / sum(w))
  q <- sum(w * (z$b_m / z$beta - b)^2)
  se <- se * max(1, sqrt(q / (nrow(z) - 1)))
  data.table(exposure = ex, mediator = med, b = b, se = se,
             p = 2 * pnorm(-abs(b / se)), n = nrow(z))
}
S1 <- rbindlist(lapply(c("GD", "AIH"), function(ex)
  rbindlist(lapply(c("c", "k"), function(md) step1(ex, md)))))

## ---- step 2: mediator -> outcome, cis instruments only --------------------
step2 <- function(oc, med) {
  o <- read9(sprintf("data/mr/cis_%s.tsv", oc))
  m <- read9(sprintf("data/mr/cis_%s.tsv", ifelse(med == "c", "CRP", "CREAT")))
  d <- merge(o[, .(key, ea, nea, bo = beta, seo = se)],
             m[, .(key, ea, nea, bm2 = beta)], by = "key", suffixes = c("", ".m"))
  d[, b_m := al(ea, nea, ea.m, nea.m, bm2)]
  d <- d[!is.na(b_m)]
  w <- 1 / d$seo^2
  b <- sum(w * d$b_m * d$bo) / sum(w * d$b_m^2)
  se <- sqrt(1 / sum(w * d$b_m^2))
  if (nrow(d) > 1) {
    q <- sum(w * (d$bo - b * d$b_m)^2)
    se <- se * max(1, sqrt(q / (nrow(d) - 1)))
  }
  data.table(outcome = oc, mediator = med, b = b, se = se, n = nrow(d))
}
OUT <- c("retic_count", "retic_pct", "Hb", "Hct", "RBC", "MCV", "MRV", "RDW", "TBIL", "DBIL")
S2 <- rbindlist(lapply(OUT, function(oc) rbindlist(lapply(c("c", "k"), function(md) step2(oc, md)))))

## ---- product of coefficients -------------------------------------------------
P <- fread("results/primary_summary.tsv")
setnames(P, c("b", "se"), c("btot", "setot"))
M <- merge(S1, S2, by = "mediator", allow.cartesian = TRUE)
M <- M[, .(exposure, outcome, mediator, b1 = b.x, se1 = se.x, b2 = b.y, se2 = se.y)]
M <- merge(M, P[, .(exposure, outcome, btot, setot)], by = c("exposure", "outcome"))
M[, ind := b1 * b2]
M[, se_ind := sqrt(b1^2 * se2^2 + b2^2 * se1^2)]
M[, prop := ind / btot]
M[, se_prop := abs(prop) * sqrt((se_ind / ind)^2 + (setot / btot)^2)]
M[, ind_lo := ind - 1.96 * se_ind][, ind_hi := ind + 1.96 * se_ind]
M[, prop_lo := prop - 1.96 * se_prop][, prop_hi := prop + 1.96 * se_prop]
setorder(M, exposure, mediator, outcome)
fwrite(S1, "results/mediation_step1.tsv", sep = "\t")
fwrite(S2, "results/mediation_step2.tsv", sep = "\t")
fwrite(M[, .(exposure, outcome, mediator, b1, se1, b2, se2, ind, se_ind, ind_lo, ind_hi,
             btot, setot, prop, prop_lo, prop_hi)],
       "results/mediation.tsv", sep = "\t")

print(S1)
print(S2[, .(outcome, mediator, b2 = round(b, 4), se = round(se, 4),
             p = signif(2 * pnorm(-abs(b / se)), 3))])
print(M[, .(exposure, outcome, mediator, indirect = round(ind * ln2, 4),
            pct = round(prop * 100, 1), lo = round(prop_lo * 100, 1),
            hi = round(prop_hi * 100, 1))])
