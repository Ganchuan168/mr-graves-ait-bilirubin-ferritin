## 02_build.R -- harmonise exposure and outcome datasets into analysis-ready form
source("00_config.R")
suppressMessages({ library(data.table); library(TwoSampleMR) })

bim <- fread(paste0(LD_REF, ".bim"), header = FALSE)
setnames(bim, c("chr19", "SNP", "cm", "pos19", "a1", "a2"))

read_expo <- function(nm) {
  d <- fread(sprintf("data/mr/instr_%s.tsv", nm))
  ## 01_instruments.R stores hg19 coordinates; older cached files are completed
  ## from the LD panel here.
  if (!all(c("chr19", "pos19") %in% names(d))) {
    d <- merge(d, bim[, .(SNP, chr19, pos19)], by = "SNP")
  }
  d[, key19 := paste(chr19, pos19, sep = ":")]
  d[, key38 := paste(chr, pos, sep = ":")]
  d
}

exposures <- list(GD = "Graves disease", AIH = "Autoimmune hypothyroidism",
                  AIT = "Autoimmune hyperthyroidism")

hg38_traits <- c("ferritin", "TSAT")

for (ex in names(exposures)) {
  e <- read_expo(ex)
  for (f in list.files("data/mr", pattern = "^out_.*tsv$", full.names = TRUE)) {
    oc <- sub("^out_", "", sub("\\.tsv$", "", basename(f)))
    ## Serum creatinine is a mediator, never a study outcome: it enters the
    ## multivariable and two-step models through out_CREAT.tsv, not through a
    ## harmonised exposure-outcome pair.
    if (oc %in% NON_OUTCOME_TRAITS) next
    o <- fread(f)
    if (nrow(o) == 0) next
    setnames(o, c("trait", "chr", "pos", "ea", "nea", "beta", "se", "eaf", "p"))
    o[, key := paste(chr, pos, sep = ":")]
    use38 <- oc %in% hg38_traits
    ekey <- if (use38) "key38" else "key19"
    m <- merge(e[, .(SNP, key = get(ekey), chr, pos, chr19, pos19,
                     beta, se, ea, nea, eaf, pval)],
               o[, .(key, beta.o = beta, se.o = se, ea.o = ea, nea.o = nea, eaf.o = eaf)],
               by = "key")
    if (nrow(m) < 3) { cat("SKIP", ex, oc, nrow(m), "\n"); next }
    d <- data.frame(SNP = m$SNP, chr = m$chr, pos = m$pos,
                    chr19 = m$chr19, pos19 = m$pos19,
                    beta.exposure = m$beta, se.exposure = m$se,
                    effect_allele.exposure = m$ea, other_allele.exposure = m$nea,
                    eaf.exposure = m$eaf, pval.exposure = m$pval,
                    beta.outcome = m$beta.o, se.outcome = m$se.o,
                    effect_allele.outcome = m$ea.o, other_allele.outcome = m$nea.o,
                    eaf.outcome = m$eaf.o,
                    exposure = exposures[[ex]], outcome = oc,
                    id.exposure = ex, id.outcome = oc, mr_keep = TRUE)
    d <- d[!is.na(d$beta.outcome) & !is.na(d$se.outcome) & d$se.outcome > 0, ]
    dir.create("data/mr/h", showWarnings = FALSE)
    fwrite(d, sprintf("data/mr/h/%s_%s.tsv", ex, oc), sep = "\t")
    cat(sprintf("%-4s %-12s merged %d SNPs\n", ex, oc, nrow(d)))
  }
}
