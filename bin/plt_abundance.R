#!/usr/bin/env Rscript

library("phyloseq")
library("tidyverse")

# parameters -------------------------------------------------------------------
# 1: phyloseq rds
# 2: significant bugs csv 
args <- commandArgs(trailingOnly = TRUE)

qiime <- readRDS(args[[1]])
sigbugs <- read.csv(args[[2]], stringsAsFactors = FALSE)

qiime_ra <- phyloseq::transform_sample_counts(qiime, function(x) x / sum(x))
phyla <- phyloseq::tax_glom(qiime_ra, "Phylum")
fam <- phyloseq::tax_glom(qiime_ra, "Family")

phyla_tax <- data.frame(tax_table(phyla)[, "Phylum"]) %>%
  tibble::rownames_to_column("taxID")
fam_tax <- data.frame(tax_table(fam)[, "Family"]) %>%
  tibble::rownames_to_column("taxID")

cohort <- data.frame(sample_data(phyla)[, "cohort"]) %>%
  tibble::rownames_to_column("sample") %>%
  mutate(sample = paste0("X", sample))

phyla_df <- data.frame(otu_table(phyla)) %>%
  tibble::rownames_to_column("taxID") %>%
  tidyr::pivot_longer(cols = X1000:X999, names_to = "sample") %>%
  left_join(phyla_tax) %>%
  left_join(cohort)

top_phyla <- phyla_df %>%
  select(value, Phylum) %>%
  group_by(Phylum) %>%
  summarise(sum = sum(value)) %>%
  arrange(desc(sum)) %>%
  top_n(9) %>%
  pull(Phylum) 

phyla_df$Phylum <- fct_other(phyla_df$Phylum, keep = top_phyla)

cohort_labeller <- function(variable,value){
  cohort_names <- list(
    'Depression'="Depression (n = 40)",
    'Healthy'="Healthy (n = 43)"
  )
  
  return(cohort_names[value])
}

ggplot(phyla_df, aes(x = sample, y = value, fill = Phylum)) +
  geom_bar(stat = "identity") + 
  facet_grid(~cohort, scales = "free", 
             labeller=cohort_labeller) + 
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  xlab("Sample") +
  ylab("Relative abundance") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(labels = scales::percent)

ggsave("phyla_abundance.png", width = 10, device = "png")

fam_df <- data.frame(otu_table(fam)) %>%
  tibble::rownames_to_column("taxID") %>%
  tidyr::pivot_longer(cols = X1000:X999, names_to = "sample") %>%
  left_join(fam_tax) %>%
  left_join(cohort)

top_family <- fam_df %>%
  select(value, Family) %>%
  group_by(Family) %>%
  summarise(sum = sum(value)) %>%
  arrange(desc(sum)) %>%
  top_n(9) %>%
  pull(Family) 

fam_df$Family <- droplevels(fct_other(fam_df$Family, keep = top_family))

ggplot(fam_df, aes(x = sample, y = value, fill = Family)) +
  geom_bar(stat = "identity") + 
  facet_grid(~cohort, scales = "free", 
             labeller=cohort_labeller) + 
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  xlab("Sample") +
  ylab("Relative abundance") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(labels = scales::percent)

ggsave("family_abundance.png", width = 10, device = "png")

sigbugs_ps <- phyloseq::prune_taxa(sigbugs$names, qiime_ra)
sigbugs_g <- data.frame(phyloseq::tax_table(sigbugs_ps)[, "Genus"]) %>%
  tibble::rownames_to_column("taxID") 
sigbugs_tax <- data.frame(phyloseq::tax_table(sigbugs_ps)[, "Species"]) %>%
  tibble::rownames_to_column("taxID") %>%
  left_join(sigbugs_g) %>%
  mutate(tax = ifelse(is.na(Species) | Species == "uncultured bacterium", paste(Genus), paste(Species))) %>%
  mutate(tax = make.unique(tax, sep = "_"))

sig_df <- data.frame(otu_table(sigbugs_ps)) %>%
  tibble::rownames_to_column("taxID") %>%
  tidyr::pivot_longer(cols = X1000:X999, names_to = "sample") %>%
  left_join(sigbugs_tax) %>%
  left_join(cohort) %>%
  filter(value != 0) %>% # long tailed data ruins box plots 
  group_by(tax) %>%
  mutate(y_max = max(value) + 0.01)

ggplot(sig_df, aes(x = cohort, y = value)) +
  geom_boxplot() +
  facet_wrap(~tax, scales = "free_y") + 
  theme_bw() +
  xlab("Sample") +
  ylab("Relative abundance") + 
  geom_jitter(position=position_jitter(0.2)) +
  xlab("") +
  ylab("Relative abundance") + 
  scale_y_continuous(labels = scales::percent) 
  # geom_blank(aes(y = y_max)) +
  # geom_signif(margin_top = 0.0001, comparisons = list(c("Depression", "Healthy"))) 

ggsave("species_abundance.png", width = 10, device = "png")