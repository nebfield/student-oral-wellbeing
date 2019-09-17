#!/usr/bin/env Rscript

library("tidyverse")

# parameters -------------------------------------------------------------------
# 1: phyloseq rds
# 2: representative sequences (FeatureData[Sequence])
args <- commandArgs(trailingOnly = TRUE)

ps <- readRDS(args[[1]])
ps_trim <- phyloseq::subset_samples(ps, smoking != "Uncertain")

rep_seqs <- qiime2R::read_qza(args[[2]])
bss_to_df <- function(dss) {
  require("Biostrings")
  return(data.frame(width=width(dss), seq=as.character(dss), names=names(dss)))
}
seq_hash <- bss_to_df(rep_seqs$data)

# check prevalence -------------------------------------------------------------
prev <- apply(X = phyloseq::otu_table(ps_trim), MARGIN =
  ifelse(phyloseq::taxa_are_rows(ps_trim), yes = 1, no = 2), FUN = function(x) {
  sum(x > 0) })

# log transform for x axis on prevalence graphs
ps_log <- phyloseq::transform_sample_counts(ps_trim, function(x) log(1 + x))

# add taxonomy and total read counts 
prevdf <- data.frame(Prevalence = prev, TotalAbundance =
  phyloseq::taxa_sums(ps_trim), phyloseq::tax_table(ps_trim))

ggplot(prevdf, aes(TotalAbundance, Prevalence / phyloseq::nsamples(ps_trim),
                    colour = Phylum)) + geom_hline(yintercept = 0.10, alpha =
                    0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
                    xlab("Total Abundance") + ylab("Prevalence") +
                    scale_y_continuous(labels = scales::percent, limits = c(0,
                    1)) + facet_wrap(~Phylum) + theme(legend.position = "none")
ggsave("prev.png", device = "png", width = 7)

# trim taxa to preserve statistical power

prevalence_threshold <- 0.10 * phyloseq::nsamples(ps_trim)
keep_taxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalence_threshold)]
ps_prev <- phyloseq::prune_taxa(keep_taxa, ps_trim)

cohortdds <- phyloseq::phyloseq_to_deseq2(ps_prev, ~ sex + smoking + cohort)

# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

geoMeans <- apply(DESeq2::counts(cohortdds), 1, gm_mean)
cohortdds <- DESeq2::estimateSizeFactors(cohortdds, geoMeans = geoMeans)
cohortdds <- DESeq2::DESeq(cohortdds, fitType="local", quiet = TRUE)
res <- DESeq2::results(cohortdds)
res <- res[order(res$padj, na.last=NA), ]
sigtab <- res[(res$padj < 0.05), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(phyloseq::tax_table(ps_prev)[rownames(sigtab), ], "matrix"))
dim(sigtab)[1] # how many bugs do we have?
sig_tidy <- sigtab %>%
  select(baseMean, log2FoldChange, padj, Phylum, Genus, Species) %>%
  tibble::rownames_to_column("names") %>% 
  left_join(seq_hash) %>%
  select(-width)
write.csv(sig_tidy, "diff_abund.csv", row.names = FALSE, quote = FALSE) 

sigtabgen <- sig_tidy
# Phylum order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$Phylum <- factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$Genus <- factor(as.character(sigtabgen$Genus), levels=names(x))
sigtabgen$Phylum <-
  factor(as.character(sigtabgen$Phylum, levels = sort(levels(sigtabgen$Phylum))))

ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) 

ggsave("diff_abund.png", device = "png", width = 10)


