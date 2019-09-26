#!/usr/bin/env Rscript

library("tidyverse")

# parameters -------------------------------------------------------------------
# 1: copyadjusted phyloseq rds
args <- commandArgs(trailingOnly = TRUE)

ps_data <- readRDS(args[[1]])

ps_data <- phyloseq::subset_samples(ps_data, smoking != "Uncertain" & smoking != "Occasional")
pslog <- phyloseq::transform_sample_counts(ps_data, function(x) log(1 + x))
# trim leftover bugs 
pslog <- phyloseq::filter_taxa(pslog, function(x) mean(x) > 0, prune = TRUE)

ord_log_cca <- phyloseq::ordinate(
  pslog,
  method = "CCA",
  distance = "bray",
  formula = pslog ~ cohort + Condition(smoking)
)

arrowmat <- vegan::scores(ord_log_cca, display = "bp")
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = CCA1, yend = CCA2, x = 0, y = 0, color = NULL)
label_map <- aes(x = 1.2 * CCA1, y = 1.2 * CCA2, color = NULL, label = labels)
arrowhead <- arrow(length = unit(0.05, "npc"))

phyloseq::plot_ordination(pslog, ord_log_cca, color = "cohort") + 
    scale_color_brewer(palette="Set1", direction = -1)
ggsave("cohort-cca.png", device = "png", width = 10)

sink("cca.txt")
vegan::anova.cca(ord_log_cca, by = "terms")
sink()

ordu = phyloseq::ordinate(pslog, "PCoA", "bray")
phyloseq::plot_ordination(pslog, ordu, color="cohort") + 
  scale_color_brewer(palette="Set1", direction = -1)
ggsave("pcoa.png", device = "png", width = 10)
