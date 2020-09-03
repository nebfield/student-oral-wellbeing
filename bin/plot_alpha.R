#!/usr/bin/env Rscript

library("tidyverse")
library("phyloseq")

# parameters -------------------------------------------------------------------
# 1: copyadjusted phyloseq rds
args <- commandArgs(trailingOnly = TRUE)

ps_data <- readRDS("ps_first.rds")

ad <- plot_richness(ps_data, x="cohort", measures=c("Observed", "InvSimpson", "Shannon")) 
df <- ad$data # just want the raw data

ggplot(df, aes(x = cohort, y = value, fill = cohort)) + 
  geom_boxplot() + 
  facet_wrap(~variable, scales = "free") +
  geom_jitter(height = 0, width = 0.1, color = "black") +
  theme_bw() + 
  xlab("") + 
  ylab("Alpha diversity measure") + 
  labs(fill = "Cohort") +
  scale_fill_manual(values=c("#66C2A5", "#FC8D62")) +
  theme(legend.position = "none")
ggsave("alpha_diversity.png", device = "png")

# Inverse Simpson is a nicer measure of alpha diversity  
# but it's not normally distributed - right tail 
df %>%
  filter(variable == "InvSimpson") %>%
  ggplot(., aes(x = value)) +
    geom_histogram()

sink("alpha.txt")
df %>%
  filter(variable == "InvSimpson") %>%
  pull(value) %>%
  shapiro.test(.)

df %>%
  filter(variable == "InvSimpson") %>%
  select(cohort, value) -> dat


wilcox.test(value ~ cohort, data = dat)
sink() 
