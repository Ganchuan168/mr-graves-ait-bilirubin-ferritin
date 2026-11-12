## 12_variant_table.R -- export the harmonised variant-level table that
## accompanies Additional file 4. One row per variant per exposure-outcome pair.
source("00_config.R")
suppressMessages({ library(data.table); library(TwoSampleMR) })

lab <- fread("results/outcome_labels.tsv")
exlab <- c(GD = "Graves' disease", AIH = "Autoimmune hypothyroidism",
           AIT = "Autoimmune hyperthyroidism")

## hg19 coordinates come from 02_build.R. Keeping them beside the hg38
## coordinates used for matching prevents the HLA flag from comparing hg38
## positions with the GRCh37 25-35 Mb boundary.
coords_for <- function(d) {
  if (all(c("chr19", "pos19") %in% names(d))) {
    return(unique(d[, .(SNP, chr, pos, chr19, pos19)]))
  }
  bim <- fread(paste0(LD_REF, ".bim"), header = FALSE,
               col.names = c("chr19", "SNP", "cm", "pos19", "a1", "a2"))
  merge(unique(d[, .(SNP, chr, pos)]), bim[, .(SNP, chr19, pos19)], by = "SNP")
}

## Clumping annotation. PLINK clumping is greedy: a variant that survives the
## r² threshold is not removed by an earlier index and therefore becomes the
## index of its own clump. Every instrument in instr_*.tsv is thus a clump lead;
## the nearest gene is carried over from the FinnGen annotation.
ann_for <- function(ex) {
  a <- fread(sprintf("data/mr/instr_%s.tsv", ex))[, .(SNP, gene)]
  a[, clump_lead := SNP]
  unique(a, by = "SNP")
}

out <- list()
for (f in list.files("data/mr/h", pattern = "\\.tsv$", full.names = TRUE)) {
  base <- sub("\\.tsv$", "", basename(f))
  ex <- sub("_.*$", "", base); oc <- sub("^[^_]+_", "", base)
  d <- fread(f)
  co <- coords_for(d)
  ed <- d[, .(SNP, beta.exposure, se.exposure, effect_allele.exposure,
              other_allele.exposure, eaf.exposure, pval.exposure,
              exposure, id.exposure)]
  od <- d[, .(SNP, beta.outcome, se.outcome, effect_allele.outcome,
              other_allele.outcome, eaf.outcome, outcome, id.outcome)]
  h <- as.data.table(suppressMessages(harmonise_data(ed, od, action = 2)))
  h <- merge(h, co, by = "SNP", all.x = TRUE)
  h <- merge(h, ann_for(ex), by = "SNP", all.x = TRUE)
  nm <- c(GD = "Graves disease", AIH = "Autoimmune hypothyroidism",
          AIT = "Autoimmune hyperthyroidism")[ex]
  n_tot <- AITD_CASES[[nm]] + AITD_CONTROLS[[nm]]
  h[, `:=`(exposure_label = exlab[ex],
           outcome_label = lab$name[match(oc, lab$outcome)],
           outcome_code = oc,
           F = (beta.exposure / se.exposure)^2,
           hla = !is.na(chr19) & chr19 == 6 & pos19 >= HLA_START & pos19 <= HLA_END,
           palindromic_flag = palindromic,
           reason = fifelse(mr_keep, "", fifelse(remove,
             "allele mismatch or missing data",
             fifelse(ambiguous, "palindromic, intermediate allele frequency",
                     "not retained"))))]
  ## F must exist before it is used: inside a single := call the symbol F would
  ## resolve to the base logical constant.
  h[, r2_approx := F / (F + n_tot - 2)]
  out[[base]] <- h[, .(exposure_label, outcome_label, outcome_code, SNP,
                       chr, pos, chr19, pos19,
                       effect_allele.exposure, other_allele.exposure, eaf.exposure,
                       beta.exposure, se.exposure, pval.exposure,
                       effect_allele.outcome, other_allele.outcome, eaf.outcome,
                       beta.outcome, se.outcome, F, r2_approx, hla,
                       palindromic_flag, mr_keep, reason, gene, clump_lead)]
}
V <- rbindlist(out)
setorder(V, exposure_label, outcome_label, SNP)
fwrite(V, "results/harmonised_variants.tsv", sep = "\t")
cat("variants:", nrow(V), " retained:", sum(V$mr_keep), "\n")
