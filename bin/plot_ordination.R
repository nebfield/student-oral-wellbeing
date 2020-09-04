#!/usr/bin/env Rscript

library("tidyverse")

# parameters -------------------------------------------------------------------
# 1: copyadjusted phyloseq rds
args <- commandArgs(trailingOnly = TRUE)

ps_data <- readRDS(args[[1]])

# bin age
phyloseq::sample_data(ps_data)$age_bin <- cut(phyloseq::sample_data(ps_data)$age, breaks = 3)

ps_data <- phyloseq::subset_samples(ps_data, smoking != "Uncertain" & smoking != "Occasional")
pslog <- phyloseq::transform_sample_counts(ps_data, function(x) log(1 + x))

# trim leftover bugs 
pslog <- phyloseq::filter_taxa(pslog, function(x) mean(x) > 0, prune = TRUE)

ord_log_cca <- phyloseq::ordinate(
  pslog,
  method = "CCA",
  distance = "bray",
  formula = pslog ~ cohort + smoking
)

ord_plot <- phyloseq::plot_ordination(pslog, ord_log_cca, color = "cohort", shape = "cohort")
arrowmat <- vegan::scores(ord_log_cca, display = "bp")
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = CCA1, yend = CCA2, x = 0, y = 0, color = NULL, shape = NULL)
label_map <- aes(x = 1.2 * CCA1, y = 1.2 * CCA2, color = NULL, shape = NULL, label = labels)
arrowhead <- arrow(length = unit(0.05, "npc"))

ord_plot + 
  scale_color_manual(values=c("#66c2a5", "#fc8d62")) + 
  geom_segment(
    arrow_map,
    size = 1.5,
    data = arrowdf,
    color = "black",
    arrow = NULL # arrowhead
  ) + 
  geom_text(label_map, size = 4, data = arrowdf) +
  stat_ellipse(type = "norm", linetype = 2) +
  labs(color='Cohort', shape = "Cohort") +
  theme_bw() 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) 
  # scale_shape_manual(values=c(19, 4)) # 19: circle, 4: cross

# Ordination ellipse https://github.com/joey711/phyloseq/issues/323

ggsave("cohort-cca.png", device = "png", dpi = 600)
ggsave("cohort-cca.svg", device = "svg")


# terms
sink("cca.txt")
vegan::anova.cca(ord_log_cca) # overall solution significant
vegan::anova.cca(ord_log_cca, by = "terms")
sink()

ordu <- phyloseq::ordinate(pslog, "PCoA", "bray")
phyloseq::plot_ordination(pslog, ordu, color = "cohort") + 
  scale_color_manual(values=c("#66C2A5", "#FC8D62")) +
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw() + 
  labs(color = "Cohort")
ggsave("pcoa.png", device = "png", dpi = 600)
ggsave("pcoa.svg", device = "svg")
