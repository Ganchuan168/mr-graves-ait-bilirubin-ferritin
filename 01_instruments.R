## 01_instruments.R -- select genetic instruments for AITD subtypes from FinnGen R12
source("00_config.R")
suppressMessages({
  library(data.table)
  library(ieugwasr)
})

bin <- PLINK_BIN
bfile <- LD_REF
refbim <- as.data.table(fread(paste0(bfile, ".bim"), header = FALSE))
setnames(refbim, c("chr", "snp", "cm", "pos", "a1", "a2"))

read_hits <- function(path, label) {
  ## column 7 of the FinnGen summary statistic is pval; pre-filter to keep the
  ## stream small, then apply the formal instrument threshold below
  d <- fread(cmd = sprintf("zcat %s | awk 'NR==1 || $7<%g'", path, P_PREFILTER),
             header = TRUE)
  setnames(d, 1:7, c("chr", "pos", "ref", "alt", "rsids", "gene", "pval"))
  d[, rsid := tstrsplit(rsids, ",", fixed = TRUE)[[1]]]
  d <- d[rsid %like% "^rs"]
  d[, `:=`(SNP = rsid, pval = as.numeric(pval), beta = as.numeric(beta),
           se = as.numeric(sebeta), ea = alt, nea = ref,
           eaf = as.numeric(af_alt), ncase = NA_real_)]
  d[, fstat := (beta / se)^2]
  d <- d[!is.na(beta) & !is.na(se) & se > 0]
  cat(label, ": SNPs with P <", P_PREFILTER, ":", nrow(d), "\n")
  ## restrict to variants present in the LD reference
  d <- d[SNP %in% refbim$snp]
  cat(label, ": present in 1KG EUR:", nrow(d), "\n")
  d
}

clump_it <- function(d, label) {
  sig <- d[pval < P_INSTRUMENT]
  cl <- ld_clump(data.frame(rsid = sig$SNP, pval = sig$pval, id = label),
                 plink_bin = bin, bfile = bfile,
                 clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2)
  out <- sig[SNP %in% cl$rsid]
  ## Every retained variant must clear the conventional weak-instrument bar.
  ## The filter is explicit so that the F-statistic counts in Table 1 and in
  ## the Results text can be reproduced from the saved instrument files.
  out <- out[fstat > 20]
  stopifnot(nrow(out) > 0, all(out$fstat > 20))
  cat(label, ": independent instruments:", nrow(out),
      "; F > 20:", sum(out$fstat > 20), "; F > 30:", sum(out$fstat > 30),
      "; F > 50:", sum(out$fstat > 50), "\n")
  out
}

files <- list(
  GD  = list(path = "data/raw/fg_GRAVES.gz",        label = "Graves disease"),
  AIH = list(path = "data/raw/fg_HYTHY_AI.gz",      label = "Autoimmune hypothyroidism"),
  AIT = list(path = "data/raw/fg_AIHYPERTHY.gz",    label = "Autoimmune hyperthyroidism")
)

res <- list()
for (nm in names(files)) {
  f <- files[[nm]]
  if (!file.exists(f$path)) { cat("missing", f$path, "\n"); next }
  d <- read_hits(f$path, f$label)
  res[[nm]] <- clump_it(d, f$label)
}

dir.create("data/mr", showWarnings = FALSE)
for (nm in names(res)) {
  x <- res[[nm]]
  ## Keep the hg19 position beside the FinnGen (hg38) position. Downstream
  ## scripts then match Pan-UK Biobank (hg19) outcomes without needing the
  ## LD panel, and the two genome builds cannot be confused.
  x <- merge(x, refbim[, .(SNP = snp, chr19 = chr, pos19 = pos)],
             by = "SNP", all.x = TRUE)
  fwrite(x[, .(SNP, chr, pos, chr19, pos19, ea, nea, eaf, beta, se, pval,
               fstat, gene)],
         sprintf("data/mr/instr_%s.tsv", nm), sep = "\t")
  cat(nm, "saved", nrow(x), "instruments; mean F =", round(mean(x$fstat), 1),
      "; min F =", round(min(x$fstat), 1), "\n")
}
