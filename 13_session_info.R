## 13_session_info.R -- record the session that actually performs the analysis.
## Running this after the pipeline loads every analysis package in the same
## process, so sessionInfo() lists TwoSampleMR, MVMR and the other dependencies
## rather than base R alone.
source("00_config.R")
suppressMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  library(MVMR)
  library(sandwich)
  library(ggplot2)
  library(cowplot)
})
sink("sessionInfo.txt")
cat("AITD-anaemia MR pipeline -- recorded after loading all analysis packages\n")

print(sessionInfo())
sink()
cat("sessionInfo.txt written\n")
