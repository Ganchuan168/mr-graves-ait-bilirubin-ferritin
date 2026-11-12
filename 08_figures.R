## 08_figures.R -- manuscript figures
source("00_config.R")
suppressMessages({ library(data.table); library(ggplot2); library(TwoSampleMR) })
dir.create("figures", showWarnings = FALSE)
ln2 <- log(2)

lab <- data.table(
  outcome = c("retic_count","retic_pct","IRF","Hb","Hct","RBC","MCV","MCH","MCHC","RDW","MRV",
         "TBIL","DBIL","ferritin","TSAT","CRP"),
  name = c("Reticulocyte count","Reticulocyte percentage","Immature reticulocyte fraction",
           "Haemoglobin","Haematocrit","Red blood cell count","Mean corpuscular volume",
           "Mean corpuscular haemoglobin","Mean corpuscular haemoglobin concentration",
           "Red cell distribution width","Mean reticulocyte volume",
           "Total bilirubin","Direct bilirubin","Ferritin","Transferrin saturation",
           "C-reactive protein"),
  domain = c(rep("Erythropoietic output",3), rep("Red cell indices",8),
             rep("Haemolysis / heme catabolism",2), rep("Iron status",2), "Inflammation"))
lab[, dom := factor(domain, levels = c("Erythropoietic output","Red cell indices","Haemolysis / heme catabolism","Iron status","Inflammation"))]

## ---- Figure 2: subtype contrast (GD minus AIH) ----
w <- fread("results/subtype_contrast.tsv")
w <- merge(w, lab, by = "outcome")
w[, `:=`(d = diff * ln2, lo = (diff - 1.96 * se_diff) * ln2, hi = (diff + 1.96 * se_diff) * ln2)]
setorder(w, dom, outcome)
w[, y := rev(seq_len(.N))]
w[, sig := p_diff < 0.05]

f2 <- ggplot(w, aes(x = d, y = y)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = "#4A5568", linewidth = 0.5) +
  geom_point(aes(fill = sig), shape = 21, size = 2.2, colour = "#4A5568", stroke = 0.5) +
  scale_fill_manual(values = c(`TRUE` = "#C53030", `FALSE` = "white"),
                    labels = c("P \u2265 0.05", "P < 0.05"), name = NULL) +
  scale_y_continuous(breaks = w$y, labels = w$name, expand = expansion(add = 0.7)) +
  labs(x = "Difference (SD per doubling; Graves minus hypothyroidism)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 8),
        legend.position = "bottom")

ggsave("figures/fig2_contrast.png", f2, width = 7.0, height = 5.2, dpi = 400)
ggsave("figures/fig2_contrast.pdf", f2, width = 6.6, height = 5.2)

## ---- Figure 3: scatter panels ----
sc <- function(ex, oc, title) {
  d <- fread(sprintf("data/mr/h/%s_%s.tsv", ex, oc))
  exd <- d[, c("SNP","beta.exposure","se.exposure","effect_allele.exposure",
               "other_allele.exposure","eaf.exposure","pval.exposure","exposure","id.exposure"), with=FALSE]
  ocd <- d[, c("SNP","beta.outcome","se.outcome","effect_allele.outcome",
               "other_allele.outcome","eaf.outcome","outcome","id.outcome"), with=FALSE]
  h <- harmonise_data(exd, ocd, action = 2); h <- h[h$mr_keep, ]
  fit <- lm(beta.outcome ~ beta.exposure, data = h, weights = 1 / h$se.outcome^2)
  cf <- coef(fit)
  ggplot(h, aes(beta.exposure, beta.outcome)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(size = 1.4, alpha = 0.7, colour = "#2B6CB0") +
    geom_abline(intercept = 0, slope = cf[2], colour = "#C53030", linewidth = 0.6) +
    labs(x = "SNP effect on exposure (log OR)", y = "SNP effect on outcome (SD)",
         title = title) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, face = "bold"))
}
g1 <- sc("GD", "retic_count", "A  Graves' disease and reticulocyte count")
g2 <- sc("AIH", "TBIL", "B  Autoimmune hypothyroidism and total bilirubin")
g <- cowplot::plot_grid(g1, g2, nrow = 1)
ggsave("figures/fig3_scatter.png", g, width = 7.6, height = 3.6, dpi = 400)
ggsave("figures/fig3_scatter.pdf", g, width = 7.6, height = 3.6)
cat("figures done\n")
