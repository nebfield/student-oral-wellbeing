#!/usr/bin/env Rscript

library("tidyverse")

# parameters
# 1) phyloseq object
# 2) picrust biom output

args <- commandArgs(trailingOnly = TRUE)

# read phyloseq data
ps_data <- readRDS(args[[1]])
ps <- phyloseq::subset_samples(ps_data, smoking != "Occasionally" & smoking != "Uncertain")

df_sd <- data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE)
lefse_meta <- df_sd %>%
  select(cohort) %>%
  tibble::rownames_to_column("subject_id") %>%
  mutate_if(is.factor, as.character) 

# read functional data
kegg <- biomformat::read_biom(args[[2]])
func <- t(as.matrix(biomformat::biom_data(kegg)))
relative_func <- prop.table(func, margin = 1) 
func_df <- as.data.frame(relative_func) %>%
  tibble::rownames_to_column("subject_id")

# lefse format
# tab separated file, rows: features, columns: samples
# first two rows: metadata, third row: subject ID

lefse <- lefse_meta %>%
  left_join(func_df, by = "subject_id") %>%
  t(.)

write.table(
  lefse,
  file = "lefse.txt",
  quote = FALSE,
  col.names = FALSE,
  row.names = TRUE,
  sep = "\t"
)
