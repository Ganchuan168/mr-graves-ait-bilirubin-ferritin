source("00_config.R")
suppressMessages({ library(data.table); library(ggplot2) })
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
p <- fread("results/primary_summary.tsv")
p <- merge(p, lab, by = "outcome")
p[, `:=`(b_d = b * ln2, lo = (b - 1.96 * se) * ln2, hi = (b + 1.96 * se) * ln2)]
p[, exposure := factor(exposure, levels = c("GD","AIH","AIT"),
                       labels = c("Graves' disease","Autoimmune hypothyroidism",
                                  "Autoimmune hyperthyroidism"))]
p[, sig := p < 0.05]
ord <- lab[, .(outcome, name, domain, dom)]; ord <- ord[outcome %in% unique(p$outcome)]
setorder(ord, dom); ord[, y := rev(seq_len(.N))]
p <- merge(p, ord[, .(outcome, name, y)], by = "outcome")
dom_col <- c("Erythropoietic output" = "#2B6CB0", "Red cell indices" = "#2F855A",
             "Haemolysis / heme catabolism" = "#B7791F", "Iron status" = "#9B2C2C",
             "Inflammation" = "#6B46C1")
p[, col := dom_col[domain]]
f1 <- ggplot(p, aes(x = b_d, y = y)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = p$col, linewidth = 0.55) +
  geom_point(aes(fill = sig), shape = 21, size = 2.4, colour = p$col, stroke = 0.5) +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"),
                    labels = c(expression(P >= 0.05), expression(P < 0.05)), name = NULL) +
  facet_wrap(~ exposure, nrow = 1) +
  scale_y_continuous(breaks = ord$y, labels = ord$name, expand = expansion(add = 0.7)) +
  labs(x = expression("Change in outcome (SD) per doubling of genetic liability"), y = NULL) +
  theme_bw(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold", size = 7.6),
        legend.position = "bottom",
        legend.text = element_text(size = 8.5),
        axis.text = element_text(size = 8.6),
        axis.text.y = element_text(size = 8.6),
        axis.title = element_text(size = 9))
ggsave("figures/fig1_forest_jtm.png", f1, width = 7.8, height = 4.9, dpi = 400)
ggsave("figures/fig1_forest_jtm.pdf", f1, width = 7.8, height = 4.9)
cat("ok\n")
