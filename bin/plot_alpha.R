#!/usr/bin/env Rscript

library("tidyverse")
library("phyloseq")

# parameters -------------------------------------------------------------------
# 1: copyadjusted phyloseq rds
args <- commandArgs(trailingOnly = TRUE)

ps_data <- readRDS(args[[1]])

ad <- plot_richness(ps_data, x="cohort", measures=c("Observed", "InvSimpson", "Shannon")) 
df <- ad$data # just want the raw data

ggplot(df, aes(x = cohort, y = value, fill = cohort)) + 
  geom_boxplot() + 
  facet_wrap(~variable, scales = "free") +
  geom_jitter(height = 0, width = 0.1, color = "black") +
  theme_bw() +
  xlab("Cohort") + 
  ylab("Alpha diversity measure") + 
  labs(fill = "Cohort") +
  scale_fill_manual(values=c("#66C2A5", "#FC8D62")) +
  theme(text = element_text(size=20, face = "bold"),
	axis.text.x = element_text(size = 12, face = "bold"))

ggsave("alpha_diversity.svg", device = "svg", width = 10, height = 5)

df %>%
  filter(variable == "Shannon") %>%
  ggplot(., aes(x = value)) +
    geom_histogram()

sink("alpha.txt")
df %>%
  filter(variable == "Shannon") %>%
  pull(value) %>%
  shapiro.test(.)

df %>%
  filter(variable == "Shannon") %>%
  select(cohort, value) -> dat

var.test(value ~ cohort, dat, alternative = "two.sided")

# Welch's t test is OK for unequal variance
t.test(value ~ cohort, data = dat)
sink() 
