## 09_mvmr3.R -- three-exposure MVMR (AITD + CRP + kidney function) and two-step mediation
source("00_config.R")
suppressMessages({ library(data.table); library(MVMR); library(sandwich) })
ln2 <- log(2)

read9 <- function(f) {
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  if (nrow(d) == 0) return(NULL)
  setnames(d, c("trait", "chr", "pos", "ea", "nea", "beta", "se", "eaf", "p"))
  d[, key := paste(chr, pos, sep = ":")]
  d
}
pull <- function(f) {
  d <- read9(f); if (is.null(d)) return(NULL)
  d[, .(key, ea, nea, beta, se)]
}
al <- function(a1, a2, ea, nea, b) {
  out <- rep(NA_real_, length(b))
  i <- ea == a1 & nea == a2
  j <- ea == a2 & nea == a1
  out[i] <- b[i]; out[j] <- -b[j]
  out
}

## ---- AITD instruments plus the cis CRP / creatinine instruments -------------
build <- function(ex, oc) {
  A <- fread(sprintf("data/mr/instr_%s.tsv", ex))[
    , .(SNP, chr19, pos19, ea, nea, beta, se)]
  A[, key := paste(chr19, pos19, sep = ":")]   # Pan-UK Biobank is hg19
  A[, `:=`(a1 = ea, a2 = nea)]
  O <- pull(sprintf("data/mr/out_%s.tsv", oc)); C <- pull("data/mr/out_CRP.tsv")
  K <- pull("data/mr/out_CREAT.tsv")
  if (is.null(O) || is.null(C) || is.null(K)) return(NULL)
  m <- merge(A, O, by = "key", suffixes = c("", ".o"))
  m <- merge(m, C, by = "key", suffixes = c("", ".c"))
  m <- merge(m, K, by = "key", suffixes = c("", ".k"))
  m[, `:=`(b_a = al(a1, a2, ea, nea, beta),         se_a = se,
           b_o = al(a1, a2, ea.o, nea.o, beta.o),   se_o = se.o,
           b_c = al(a1, a2, ea.c, nea.c, beta.c),   se_c = se.c,
           b_k = al(a1, a2, ea.k, nea.k, beta.k),   se_k = se.k)]
  mm <- m[!is.na(b_a) & !is.na(b_o) & !is.na(b_c) & !is.na(b_k),
          .(key, b_a, se_a, b_c, se_c, b_k, se_k, b_o, se_o)]

  ## cis block: 3 CRP-locus variants + 1 UMOD variant, effects on every trait
  cz <- merge(pull(sprintf("data/mr/cis_%s.tsv", oc)), pull("data/mr/cis_CRP.tsv"),
              by = "key", suffixes = c("", ".c"))
  cz <- merge(cz, pull("data/mr/cis_CREAT.tsv"), by = "key", suffixes = c("", ".k"))
  ca <- fread("data/mr/cis_aitd.tsv")[exposure == ex, .(key, ea, nea, beta, se)]
  cz <- merge(cz, ca, by = "key", suffixes = c("", ".a"))
  cz[, `:=`(b_a = al(ea, nea, ea.a, nea.a, beta.a), se_a = se.a,
            b_o = beta, se_o = se, b_c = beta.c, se_c = se.c,
            b_k = beta.k, se_k = se.k)]
  cz <- cz[!is.na(b_a) & !is.na(b_o) & !is.na(b_c) & !is.na(b_k),
           .(key, b_a, se_a, b_c, se_c, b_k, se_k, b_o, se_o)]
  list(z = rbind(mm[!key %in% cz$key], cz), ncis = nrow(cz))
}

fit_hc0 <- function(z, cols) {
  f <- as.formula(paste("b_o ~", paste(cols, collapse = " + "), "- 1"))
  fit <- lm(f, data = z, weights = 1 / z$se_o^2)
  V <- vcovCL(fit, type = "HC0"); b <- coef(fit); se <- sqrt(diag(V))
  list(b = b, se = se, p = 2 * pnorm(-abs(b / se)), n = nrow(z))
}

OUT <- c("retic_count", "retic_pct", "Hb", "Hct", "RBC", "MCV", "MRV", "RDW",
         "TBIL", "DBIL", "CRP")
res <- list()
for (ex in c("GD", "AIH")) for (oc in OUT) {
  if (oc == "CRP") next
  bz <- build(ex, oc)
  if (is.null(bz)) next
  z <- bz$z; ncis <- bz$ncis
  if (nrow(z) < 8) next
  r  <- fit_hc0(z, c("b_a", "b_c", "b_k"))
  r2 <- fit_hc0(z, c("b_a", "b_c"))
  u  <- fit_hc0(z, c("b_a"))
  inp <- format_mvmr(BXGs = as.matrix(z[, .(b_a, b_c, b_k)]), BYG = z$b_o,
                     seBXGs = as.matrix(z[, .(se_a, se_c, se_k)]), seBYG = z$se_o,
                     RSID = z$key)
  cf <- try(strength_mvmr(inp, gencov = 0), silent = TRUE)
  cfv <- if (inherits(cf, "try-error")) rep(NA_real_, 3) else as.numeric(cf[1, ])
  res[[paste(ex, oc)]] <- data.table(
    exposure = ex, outcome = oc, n = r$n, n_cis = ncis,
    bAITD = r$b["b_a"], seAITD = r$se["b_a"], pAITD = r$p["b_a"],
    bCRP = r$b["b_c"], seCRP = r$se["b_c"], pCRP = r$p["b_c"],
    bKID = r$b["b_k"], seKID = r$se["b_k"], pKID = r$p["b_k"],
    bAITD_crp = r2$b["b_a"], seAITD_crp = r2$se["b_a"], pAITD_crp = r2$p["b_a"],
    bUniv = u$b["b_a"], seUniv = u$se["b_a"], pUniv = u$p["b_a"],
    condF_AITD = cfv[1], condF_CRP = cfv[2], condF_KID = cfv[3])
}
R <- rbindlist(res)
if (nrow(R) == 0) stop("no MVMR model could be built; check data/mr/instr_*.tsv, out_*.tsv and cis_aitd.tsv")

## ---- two-step mediation ----------------------------------------------------
## step 1: AITD -> mediator, using only the AITD instruments
## step 1: IVW using the AITD instruments only (Wald ratio weighted by b_a^2/se_m^2)
step1 <- function(ex, med) {
  ## FinnGen (instr) is hg38; Pan-UK Biobank (out_*) is hg19. Use the hg19
  ## coordinates stored by 01_instruments.R, exactly as build() does.
  A <- fread(sprintf("data/mr/instr_%s.tsv", ex))[
    , .(SNP, chr19, pos19, ea, nea, beta, se)]
  A[, key := paste(chr19, pos19, sep = ":")]
  M <- fread(sprintf("data/mr/out_%s.tsv", ifelse(med == "c", "CRP", "CREAT")))
  setnames(M, c("trait", "chr", "pos", "ea", "nea", "beta", "se", "eaf", "p"))
  M[, key := paste(chr, pos, sep = ":")]
  z0 <- merge(A, M[, .(key, ea.m = ea, nea.m = nea, bm = beta, sem = se)], by = "key")
  z0[, b_m := al(ea, nea, ea.m, nea.m, bm)]
  z0 <- z0[!is.na(b_m)]
  if (nrow(z0) == 0) return(NULL)
  w <- z0$beta^2 / z0$sem^2
  if (!is.finite(sum(w)) || sum(w) <= 0) return(NULL)
  b <- sum(w * (z0$b_m / z0$beta)) / sum(w)
  se <- sqrt(1 / sum(w))
  ## Same multiplicative random-effects inflation as 10_mediation.R, so the
  ## two scripts cannot disagree about the first-step standard error.
  if (nrow(z0) > 1) {
    q <- sum(w * (z0$b_m / z0$beta - b)^2)
    se <- se * max(1, sqrt(q / (nrow(z0) - 1)))
  }
  data.table(exposure = ex, mediator = med, b = b, se = se,
             p = 2 * pnorm(-abs(b / se)), n = nrow(z0))
}
S1 <- rbindlist(lapply(c("GD", "AIH"), function(ex)
  rbindlist(lapply(c("c", "k"), function(md) step1(ex, md)))))
if (nrow(S1) == 0) stop("first-step mediation produced no rows; check hg19 keys in instr_*.tsv")

## step 2: mediator -> outcome, using only the cis instruments for the mediator
step2 <- function(oc, med) {
  o <- pull(sprintf("data/mr/cis_%s.tsv", oc)); m <- pull(sprintf("data/mr/cis_%s.tsv",
    ifelse(med == "c", "CRP", "CREAT")))
  if (is.null(o) || is.null(m)) return(NULL)
  d <- merge(o, m, by = "key", suffixes = c("", ".m"))
  d[, b_m := al(ea, nea, ea.m, nea.m, beta.m)]
  d <- d[!is.na(b_m) & !is.na(beta)]
  if (nrow(d) == 0) return(NULL)
  w <- 1 / d$se^2
  b <- sum(w * d$b_m * d$beta) / sum(w * d$b_m^2)
  se <- sqrt(1 / sum(w * d$b_m^2))
  ## multiplicative random-effects SE if >1 variant
  if (nrow(d) > 1) {
    q <- sum(w * (d$beta - b * d$b_m)^2)
    se <- se * max(1, sqrt(q / (nrow(d) - 1)))
  }
  data.table(outcome = oc, mediator = med, b = b, se = se, n = nrow(d))
}
S2 <- rbindlist(lapply(OUT, function(oc) rbindlist(lapply(c("c", "k"), function(md) {
  r <- try(step2(oc, md), silent = TRUE); if (!inherits(r, "try-error")) r
}))))

fwrite(R,  "results/mvmr3.tsv", sep = "\t")
## Distinct names: 10_mediation.R owns mediation_step1/2.tsv, and 09 must not
## overwrite them with a differently standardised first step.
fwrite(S1, "results/mvmr3_step1.tsv", sep = "\t")
fwrite(S2, "results/mvmr3_step2.tsv", sep = "\t")

print(R[, .(exposure, outcome, n, n_cis,
            univ = round(bUniv * ln2, 4), adjCRP = round(bAITD_crp * ln2, 4),
            adjBoth = round(bAITD * ln2, 4), pBoth = signif(pAITD, 3),
            F = round(condF_AITD, 1), Fc = round(condF_CRP, 1),
            Fk = round(condF_KID, 1))])
print(S1)
print(S2)
