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
trim_bugs <- function(ps, level) {
  prev <- apply(X = phyloseq::otu_table(ps), MARGIN =
                  ifelse(phyloseq::taxa_are_rows(ps), yes = 1, no = 2), FUN = function(x) {
                    sum(x > 0) })
  
  prevdf <- data.frame(
    Prevalence = prev,
    TotalAbundance =
      phyloseq::taxa_sums(ps),
    phyloseq::tax_table(ps)
  )
  
  # keep if more than 5% prevalence across all samples
  prevalence_threshold <- 0.05 * phyloseq::nsamples(ps)
  keep_taxa_prev <- rownames(prevdf)[(prevdf$Prevalence >= prevalence_threshold)]
  ps_prev <- phyloseq::prune_taxa(keep_taxa_prev, ps)
  
  # keep if more than 0.1% relative abundance 
  ps_glom <- phyloseq::tax_glom(ps_prev, level)
  keep_taxa <- phyloseq::taxa_sums(ps_glom) / sum(phyloseq::taxa_sums(ps_glom)) * 100 > 0.001
  ps_trim <- phyloseq::prune_taxa(keep_taxa, ps_glom)
  
  # now transform to RA
  return(phyloseq::transform_sample_counts(ps_trim, function(x) x / sum(x)))
}
phyla <- trim_bugs(qiime, "Phylum")
fam <- trim_bugs(qiime, "Family")

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

cohort_labeller <- function(variable,value){
  cohort_names <- list(
    'Healthy'="Healthy (n = 43)",
    'Depression'="Depression (n = 40)"
  )
  
  return(cohort_names[value])
}

phyla_df %>%
  group_by(Phylum) %>%
  summarise(total = sum(value)) %>%
  arrange(desc(total)) %>%
  pull(Phylum) -> phyla_order

# rearrange factor level
phyla_df$Phylum <- forcats::fct_relevel(phyla_df$Phylum, as.character(phyla_order))

phyla_df %>%
  group_by(cohort) %>%
  mutate(group_no = match(sample, unique(sample))) ->
  phyla_df 

ggplot(phyla_df, aes(x = group_no, y = value, fill = Phylum)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_wrap(~cohort, scales = "free",
             labeller=cohort_labeller) + 
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  xlab("Sample") +
  ylab("Relative abundance") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 44, by = 1))

ggsave("phyla_abundance.png", width = 10, device = "png")

fam_df <- data.frame(otu_table(fam)) %>%
  tibble::rownames_to_column("taxID") %>%
  tidyr::pivot_longer(cols = X1000:X999, names_to = "sample") %>%
  left_join(fam_tax) %>%
  left_join(cohort)

# top 10 
top_family <- fam_df %>%
  select(value, Family) %>%
  group_by(Family) %>%
  summarise(sum = sum(value)) %>%
  arrange(desc(sum)) %>%
  top_n(9) %>%
  pull(Family) 

fam_df$Family <- droplevels(fct_other(fam_df$Family, keep = top_family, other_level = "Other (n=40)"))

fam_df %>%
  group_by(Family) %>%
  summarise(total = sum(value)) %>%
  arrange(desc(total)) %>%
  pull(Family) -> fam_order

# rearrange factor level
fam_df$Family <- forcats::fct_relevel(fam_df$Family, as.character(fam_order))

fam_df %>%
  group_by(cohort) %>%
  mutate(group_no = match(sample, unique(sample))) ->
  fam_df 

ggplot(fam_df, aes(x = group_no, y = value, fill = Family)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_wrap(~cohort, scales = "free",
             labeller=cohort_labeller) + 
  theme_bw() +
  xlab("Sample") +
  ylab("Relative abundance") + 
  scale_y_continuous(labels = scales::percent) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 44, by = 1))

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