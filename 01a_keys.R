## 01a_keys.R -- variant position lists used by the download scripts
source("00_config.R")
suppressMessages(library(data.table))

bim <- fread(paste0(LD_REF, ".bim"), header = FALSE,
             col.names = c("chr", "SNP", "cm", "pos", "a1", "a2"))

d <- rbindlist(lapply(c("GD", "AIH", "AIT"), function(x)
  fread(sprintf("data/mr/instr_%s.tsv", x))))

## hg38: coordinates as reported by FinnGen
writeLines(sort(unique(paste(d$chr, d$pos, sep = ":"))), "data/mr/keys_hg38.txt")

## hg19: coordinates from the LD reference panel
m <- merge(d[, .(SNP)], bim[, .(SNP, chr, pos)], by = "SNP")
writeLines(sort(unique(paste(m$chr, m$pos, sep = ":"))), "data/mr/keys_hg19.txt")

cat("hg38 keys:", nrow(unique(d[, .(chr, pos)])),
    " hg19 keys:", nrow(unique(m[, .(chr, pos)])), "\n")
