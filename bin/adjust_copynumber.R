#!/usr/bin/env Rscript

library("tidyverse")

# transform abundance data to incorporate 16S copy number data
# "Incorporating 16S Gene Copy Number Information Improves Estimates of Microbial
# Diversity and Abundance" Kembel et al. 2012
# https://doi.org/10.1371/journal.pcbi.1002743

# parameters -------------------------------------------------------------------
# 1: path to phyloseq object (rdata file) ps_data
args <- commandArgs(trailingOnly = TRUE)

ps_data <- readRDS(args[[1]])

rrnDB_url <- "https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.5_pantaxa_stats_RDP.tsv.zip"
tf <- tempfile()
td <- tempdir()
download.file(rrnDB_url, tf)
unzip(tf, exdir = td)

rrnDB5_5 <- read.table(file.path(td, "rrnDB-5.5_pantaxa_stats_RDP.tsv"), sep =
  "\t", header = TRUE, stringsAsFactors = FALSE )

# get operon number of bugs present in database --------------------------------

operons <- rrnDB5_5 %>% select(name, median) %>% rename(Genus = name)
genera <- data.frame(phyloseq::tax_table(ps_data)[, "Genus"], stringsAsFactors = FALSE)
genera <- genera %>% tibble::rownames_to_column("sequence") # preserve rownames
known_operons <- genera %>%
  inner_join(operons, by = "Genus") %>%
  tibble::column_to_rownames("sequence") %>%
  select(-Genus) %>%
  picante::df2vec()

# estimate operon number for unknown bugs --------------------------------------

# How can I resolve polytomies in my phylogeny?
# https://en.wikipedia.org/wiki/Polytomy
# https://www.r-phylo.org/wiki/HowTo/DataTreeManipulation

dichotomousphylogeny <- ape::multi2di(phyloseq::phy_tree(ps_data))
estimated_operons <- picante::phyEstimate(dichotomousphylogeny, known_operons)
estimated_operons$estimate <- round(estimated_operons$estimate, digits = 1)
known_operons_df <- data.frame(copynum = known_operons)
estimated_operons_df <- data.frame(estimated_operons) %>%
  rename(copynum = estimate) %>%
  select(-se)
# combine known and estimated operon numbers
operons <- rbind(known_operons_df, estimated_operons_df)

# transform abundances and save to rdata file ----------------------------------

adj_copynum <- function(observed, copynum) {
  return(round(observed / copynum, digits = 0))
}

otu <- bind_cols(data.frame(phyloseq::otu_table(ps_data)), operons)
otu %>%
  group_by(copynum) %>%
  mutate_at(vars(-group_cols()), list(~ adj_copynum(., copynum))) %>%
  ungroup() %>%
  select(-copynum) -> adj_otu

adj_otu_ps <- phyloseq::otu_table(data.matrix(adj_otu), taxa_are_rows = TRUE)
rownames(adj_otu_ps) <- phyloseq::taxa_names(ps_data)
colnames(adj_otu_ps) <- phyloseq::sample_names(ps_data) 
ps_copyadj <- ps_data 
phyloseq::otu_table(ps_copyadj) <- phyloseq::otu_table(adj_otu_ps)

data.frame(
  orig_count = phyloseq::taxa_sums(ps_data),
  adj_count = phyloseq::taxa_sums(ps_copyadj),
  copynum = operons$copynum
) %>%
  tibble::rownames_to_column("id") %>%
  tibble::as_tibble(.) %>%
  arrange(desc(copynum)) -> adj_log

write.csv(adj_log, file = "adj_log.csv", quote = FALSE, row.names = FALSE)

# drop any taxa with all 0 reads after copy number adjustment
ps_copyadj <-
  phyloseq::prune_taxa(phyloseq::taxa_sums(ps_copyadj) > 0, ps_copyadj)
saveRDS(ps_copyadj, file = args[[2]])