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

# this function is probably far too complicated but I'm not smart enough to
# think of a simpler way to do it
adjust_copy_number <- function(x, ...) {
  args <- list(...)
  operons <- data.frame(args$operons) %>% tibble::rownames_to_column("seq")
  sample_df <- data.frame(abund = x) %>% tibble::rownames_to_column("seq")
  joined <- dplyr::left_join(sample_df, operons, by = "seq")

  # phyloseq uses dummy data to test validity of a user-defined function: 1:10
  # the join fails to add a copynum column to this dummy data, which will
  # cause the mutate below to break
  if (!("copynum" %in% names(joined))) {
    joined$copynum <- 1:10
  }

  # round to 0 digits to maintain counts
  copynum_adj <-
    dplyr::mutate(joined, abund_adj = round(abund / copynum, digits = 0)) %>%
    select(seq, abund_adj) %>%
    tibble::column_to_rownames("seq") %>%
    picante::df2vec()

  return(copynum_adj)
}

ps_copyadj <-
  phyloseq::transform_sample_counts(ps_data,
                                    adjust_copy_number,
                                    operons = operons)

# drop any taxa with all 0 reads after copy number adjustment
ps_copyadj <-
  phyloseq::prune_taxa(phyloseq::taxa_sums(ps_copyadj) > 0, ps_copyadj)
saveRDS(ps_copyadj, file = "ps_copyadj.rds")

write.table(
  operons %>% tibble::rownames_to_column("seq") %>% select(copynum, seq),
  file = "copynumber.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

