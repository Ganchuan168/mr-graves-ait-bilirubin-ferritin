## 10b_mediation_bootstrap.R -- parametric bootstrap for mediated proportions
##
## The delta-method interval for a ratio treats the indirect and total effects
## as independent and can extend outside [0, 1]. This script draws B values of
## b1, b2 and b_total from their estimated sampling distributions and reports
## percentile intervals. The draws are approximate: they ignore the covariance
## between the two steps and the uncertainty of the second-step instruments.

source("00_config.R")
suppressMessages(library(data.table))
set.seed(SEED)

m <- fread("results/mediation.tsv")
B <- N_BOOT
out <- m[, {
  b1d <- rnorm(B, b1, se1)
  b2d <- rnorm(B, b2, se2)
  btd <- rnorm(B, btot, setot)
  p <- (b1d * b2d) / btd
  list(prop_boot_lo = quantile(p, 0.025, na.rm = TRUE),
       prop_boot_hi = quantile(p, 0.975, na.rm = TRUE),
       prop_boot_median = median(p, na.rm = TRUE),
       n_boot = sum(is.finite(p)))
}, by = .(exposure, outcome, mediator)]

fwrite(out, "results/mediation_bootstrap.tsv", sep = "\t")
print(out, nrows = 100)
