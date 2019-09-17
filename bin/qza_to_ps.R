#!/usr/bin/env Rscript

library("tidyverse")

# parameters -------------------------------------------------------------------
# 1: feature table path
# 2: rooted tree path
# 3: taxonomy path
# 4: sample metadata path

args <- commandArgs(trailingOnly = TRUE)

df <- tibble::as_tibble(read.table(args[[4]], header= TRUE, sep = "\t"))

# edit sample data: replace numeric data with character factors, set levels ----
# Smoking -> 3: Occasionally, 1: Past -> Uncertain
samp_dat <- df %>%
  mutate(sex = factor(sex)) %>%
  mutate(cohort = factor(cohort)) %>%
  mutate(SMOKING = factor(SMOKING)) %>%
  mutate(id = as.character(id)) %>%
  mutate(sex = forcats::fct_recode(sex, "Male" = "1", "Female" = "2")) %>%
  mutate(smoking = forcats::fct_recode(SMOKING, "Never" = "4", "Uncertain" = "3", "Daily" = "2", "Uncertain" = "1")) %>%
  mutate(cohort = forcats::fct_recode(cohort, "Healthy" = "0", "Depression" = "1")) %>%
  mutate(cohort = relevel(cohort, "Healthy")) %>%
  select(-SMOKING, -dep_score)

write.table(samp_dat, "fsd.tsv", row.names = FALSE, sep = "\t", quote = FALSE)

ps <- qiime2R::qza_to_phyloseq(args[[1]], args[[2]], args[[3]], "fsd.tsv")

# fix tax table
# see: https://github.com/jbisanz/qiime2R/issues/14
tax <- data.frame(phyloseq::tax_table(ps)[, 1]) %>%
  mutate(Kingdom = stringr::str_replace_all(Kingdom, "D_\\d__", ""))
tax <- tax %>%
  tidyr::separate(Kingdom, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

tax_mat <- as.matrix(tax)
rownames(tax_mat) <- phyloseq::taxa_names(ps)
phyloseq::tax_table(ps) <- tax_mat

# clean samples ----------------------------------------------------------------
# problems with the following samples:
remove_samps <- !phyloseq::sample_names(ps) %in% c("238", "907", "708", "763")
ps_pruned <- phyloseq::prune_samples(remove_samps, ps)

# also remove subjects with no smoking data recorded
missing_smokers <- !is.na(phyloseq::sample_data(ps_pruned)$smoking)
ps_pruned <- phyloseq::prune_samples(missing_smokers, ps_pruned)

# remove any bugs that occur 0 times across all samples
# this could happen because some bugs unique to the removed samples might still
# be present in the tax table and phylogenetic tree of the phyloseq object
ps_pruned <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_pruned) > 0, ps_pruned)

# remove sequence variants with ambiguous phylum annotation 
# these sequence variants are probably just artefacts
ps_pruned <- phyloseq::subset_taxa(ps_pruned, !is.na(Phylum)) 

# relevel
phyloseq::sample_data(ps_pruned)$cohort <- relevel(phyloseq::sample_data(ps_pruned)$cohort, "Healthy")

saveRDS(ps_pruned, file = "ps_first.rds")
