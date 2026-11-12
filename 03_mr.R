## 03_mr.R -- two-sample MR using harmonised datasets in data/mr/h
source("00_config.R")
suppressMessages({ library(data.table); library(TwoSampleMR) })
set.seed(SEED)
dir.create("results", showWarnings = FALSE)

label <- c(
  retic_count = "Reticulocyte count", retic_pct = "Reticulocyte percentage",
  IRF = "Immature reticulocyte fraction", MRV = "Mean reticulocyte volume",
  RDW = "Red cell distribution width", MCV = "Mean corpuscular volume",
  MCH = "Mean corpuscular haemoglobin", MCHC = "MCHC",
  Hb = "Haemoglobin", RBC = "Red blood cell count", Hct = "Haematocrit",
  TBIL = "Total bilirubin", DBIL = "Direct bilirubin",
  ferritin = "Ferritin", TSAT = "Transferrin saturation", CRP = "C-reactive protein"
)
exlab <- c(GD = "Graves disease", AIH = "Autoimmune hypothyroidism",
           AIT = "Autoimmune hyperthyroidism")

files <- list.files("data/mr/h", pattern = "\\.tsv$", full.names = TRUE)
res <- list(); het <- list(); plei <- list(); loo_s <- list(); steig <- list()
for (f in files) {
  d <- fread(f)
  if (nrow(d) < 5) next
  ex <- d$id.exposure[1]; oc <- d$id.outcome[1]
  ex_dat <- d[, c("SNP","beta.exposure","se.exposure","effect_allele.exposure",
                  "other_allele.exposure","eaf.exposure","pval.exposure",
                  "exposure","id.exposure"), with = FALSE]
  oc_dat <- d[, c("SNP","beta.outcome","se.outcome","effect_allele.outcome",
                  "other_allele.outcome","eaf.outcome","outcome","id.outcome"), with = FALSE]
  h <- harmonise_data(ex_dat, oc_dat, action = 2)
  h <- h[h$mr_keep, ]
  if (nrow(h) < 5) { cat("skip", ex, oc, nrow(h), "\n"); next }
  r <- mr(h, method_list = c("mr_ivw_mre", "mr_egger_regression",
                             "mr_weighted_median", "mr_weighted_mode"))
  r$exposure <- exlab[[ex]]; r$outcome <- label[[oc]]; r$n_snp <- nrow(h)
  res[[paste(ex, oc)]] <- r

  hh <- mr_heterogeneity(h, method_list = c("mr_ivw_mre", "mr_egger_regression"))
  hh$exposure <- ex; hh$outcome <- oc; het[[paste(ex, oc)]] <- hh

  pl <- mr_pleiotropy_test(h); pl$exposure <- ex; pl$outcome <- oc
  plei[[paste(ex, oc)]] <- pl

  lo <- mr_leaveoneout(h)
  loo_s[[paste(ex, oc)]] <- data.table(exposure = ex, outcome = oc,
    loo_min = min(lo$b), loo_max = max(lo$b),
    loo_sign_stable = all(sign(lo$b) == sign(lo$b[1])))

  st <- tryCatch(directionality_test(h), error = function(e) NULL)
  if (!is.null(st)) { st$exposure <- ex; st$outcome <- oc; steig[[paste(ex, oc)]] <- st }

  ivw <- r[r$method == "Inverse variance weighted (multiplicative random effects)", ]
  cat(sprintf("%-4s %-12s n=%3d b=%+.4f se=%.4f p=%.3g\n", ex, oc, nrow(h), ivw$b, ivw$se, ivw$pval))
}
R <- rbindlist(res, fill = TRUE); fwrite(R, "results/mr_estimates.tsv", sep = "\t")
if (length(het)) fwrite(rbindlist(het, fill = TRUE), "results/heterogeneity.tsv", sep = "\t")
if (length(plei)) fwrite(rbindlist(plei, fill = TRUE), "results/pleiotropy.tsv", sep = "\t")
if (length(loo_s)) fwrite(rbindlist(loo_s, fill = TRUE), "results/leaveoneout.tsv", sep = "\t")
if (length(steig)) fwrite(rbindlist(steig, fill = TRUE), "results/steiger.tsv", sep = "\t")
cat("done\n")
