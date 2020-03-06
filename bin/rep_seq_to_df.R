#!/usr/bin/env Rscript

# parameters -------------------------------------------------------------------
# 1: representative sequences (FeatureData[Sequence])
args <- commandArgs(trailingOnly = TRUE)

rep_seqs <- qiime2R::read_qza(args[[1]])

bss_to_df <- function(dss) {
  require("Biostrings")
  return(data.frame(width=width(dss), seq=as.character(dss), names=names(dss)))
}

seq_hash <- bss_to_df(rep_seqs$data)
saveRDS(seq_hash, "seq_hash.rds")