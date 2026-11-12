## 05_mvmr.R -- multivariable MR: AITD adjusted for CRP
source("00_config.R")
suppressMessages({ library(data.table); library(TwoSampleMR); library(MVMR); library(sandwich) })

## Build MVMR input from the merged exposure-outcome files. Each pair is
## harmonised first, so that exposure and outcome effects are on the same
## allele before the two instrument blocks are combined.
harm_pair <- function(f) {
  d <- fread(f)
  ed <- d[, .(SNP, beta.exposure, se.exposure, effect_allele.exposure,
              other_allele.exposure, eaf.exposure, pval.exposure,
              exposure, id.exposure)]
  od <- d[, .(SNP, beta.outcome, se.outcome, effect_allele.outcome,
              other_allele.outcome, eaf.outcome, outcome, id.outcome)]
  h <- as.data.table(suppressMessages(harmonise_data(ed, od, action = 2)))
  h[h$mr_keep == TRUE, ]
}

build <- function(ex, oc) {
  f1 <- sprintf("data/mr/h/%s_%s.tsv", ex, oc)
  f2 <- sprintf("data/mr/h/%s_CRP.tsv", ex)
  if (!file.exists(f1) || !file.exists(f2)) return(NULL)
  a <- harm_pair(f1); b <- harm_pair(f2)
  if (nrow(a) < 5 || nrow(b) < 5) return(NULL)
  m <- merge(
    a[, .(SNP, beta.exp1 = beta.exposure, se.exp1 = se.exposure,
          beta.out = beta.outcome, se.out = se.outcome)],
    b[, .(SNP, beta.exp2 = beta.outcome, se.exp2 = se.outcome)],
    by = "SNP")
  m[!is.na(beta.out) & !is.na(beta.exp2)]
}

mvmr_ivw <- function(m) {
  w <- 1 / m$se.out^2
  fit <- lm(beta.out ~ beta.exp1 + beta.exp2 - 1, data = m, weights = w)
  V <- vcovCL(fit, type = "HC0")
  b <- coef(fit); se <- sqrt(diag(V))
  data.table(term = names(b), b = as.numeric(b), se = as.numeric(se),
             p = 2 * pnorm(-abs(as.numeric(b) / as.numeric(se))))
}

## conditional F statistics of Sanderson et al. via MVMR::strength_mvmr
cond_f <- function(m) {
  inp <- format_mvmr(
    BXGs = as.matrix(m[, .(beta.exp1, beta.exp2)]),
    BYG = m$beta.out,
    seBXGs = as.matrix(m[, .(se.exp1, se.exp2)]),
    seBYG = m$se.out,
    RSID = m$SNP)
  cf <- try(strength_mvmr(inp, gencov = 0), silent = TRUE)
  if (inherits(cf, "try-error")) return(c(NA_real_, NA_real_))
  as.numeric(cf[1, ])
}

## Same ten outcomes as the three-exposure model in 09_mvmr3.R, so the
## CRP-adjusted and fully adjusted results cover an identical outcome set.
OUT <- c("retic_count", "retic_pct", "Hb", "Hct", "RBC", "MCV", "MRV", "RDW",
         "TBIL", "DBIL")

res <- list()
for (ex in c("GD", "AIH")) {
  for (oc in OUT) {
    m <- build(ex, oc); if (is.null(m) || nrow(m) < 10) next
    r <- mvmr_ivw(m); r[, `:=`(exposure = ex, outcome = oc, n = nrow(m))]
    cf <- cond_f(m)
    r[, condF_AITD := cf[1]][, condF_CRP := cf[2]]
    res[[paste(ex, oc)]] <- r
  }
}
R <- rbindlist(res)
fwrite(R, "results/mvmr.tsv", sep = "\t")
print(R[term == "beta.exp1", .(exposure, outcome, n,
      b_AITD_adj = round(b, 4), se = round(se, 4), p = signif(p, 3),
      condF_AITD = round(condF_AITD, 1), condF_CRP = round(condF_CRP, 1))])
