#!/usr/bin/env Rscript

library("tidyverse")
library("phyloseq")

# parameters -------------------------------------------------------------------
# 1: copyadjusted phyloseq rds
args <- commandArgs(trailingOnly = TRUE)

ps_data <- readRDS(args[[1]])

ad <- plot_richness(ps_data, x="cohort", measures=c("Observed", "InvSimpson", "Shannon")) 
ad <- ad + geom_violin() +
  geom_jitter(height = 0, width = 0.1, color = "black") +
  theme_classic() + 
  xlab("")

ggplot(ad$data, aes(x = cohort, y = value, fill = cohort)) + 
  geom_violin() + 
  facet_wrap(~variable, scales = "free") +
  geom_jitter(height = 0, width = 0.1, color = "black") +
  theme_linedraw() + 
  xlab("") + 
  ylab("Alpha diversity measure") + 
  labs(fill = "Cohort") +
  scale_fill_manual(values=c("#66C2A5", "#FC8D62")) 
  # scale_fill_brewer(palette="Set2", direction = 1) 
  # scale_fill_manual(values=c("#B2DF8A", "#FB9A99"))

ggsave("alpha_diversity.png", device = "png", width = 10)

d <- ad$data %>%
  group_split(variable, cohort) 

sink("alpha.txt")
print(d)
t.test(d[[1]]$value, d[[2]]$value)
t.test(d[[3]]$value, d[[4]]$value)
t.test(d[[5]]$value, d[[6]]$value)
sink()
