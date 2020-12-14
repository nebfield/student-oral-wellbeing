#!/usr/bin/env Rscript

library("tidyverse")

# parameters -------------------------------------------------------------------
# 1: phyloseq rds
# 2: representative sequences dataframe
args <- commandArgs(trailingOnly = TRUE)

ps <- readRDS(args[[1]])
ps_trim <- phyloseq::subset_samples(ps, smoking != "Uncertain")
seq_hash <- readRDS(args[[2]])

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

sig_tidy %>%
  mutate(Species = na_if(Species, "uncultured bacterium")) %>%
  mutate(lab = ifelse(is.na(Species), paste0(Genus), paste0(Species))) %>%
  mutate(lab = make.unique(lab, sep = " ")) %>%
  mutate(Genus = as.factor(Genus)) %>%
  tibble::as_tibble(.) -> sigtabgen

# colour palette for loads of data
# https://www.r-bloggers.com/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/
# colourCount = length(unique(sigtabgen$Genus))
# getPalette = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))

# order bars from lowest to highest
sigtabgen %>%
  arrange(desc(log2FoldChange)) %>%
  pull(lab) -> sig_order

# rearrange x labels 
sigtabgen$lab <- forcats::fct_relevel(sigtabgen$lab, as.character(sig_order))

# rearrange legend form alphabetical to order of appearance (phylum) 
sigtabgen$Phylum <-
  forcats::fct_relevel(
    sigtabgen$Phylum,
    c(
      "Proteobacteria",
      "Actinobacteria",
      "Bacteroidetes",
      "Spirochaetes",
      "Firmicutes",
      "Fusobacteria"
    ))
    
ggplot(sigtabgen, aes(y=log2FoldChange, x=lab, fill=Phylum)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
  coord_flip() +
  theme_bw() + 
  xlab("") + 
  ylab(expression(Log["2"]*"-fold change")) + 
  # scale_fill_brewer(palette = "Paired") +
  scale_fill_manual(values=c("#B2DF8A", "#33A02C", "#A6CEE3", "#FF7F00", "#1F78B4", "#FB9A99")) +  
  theme(text = element_text(size=12, face ="bold"),
        axis.text.y = element_text(face = "italic"))
  
ggsave("diff_abund.png", device = "png", width = 10)
ggsave("diff_abund.svg", device = "svg", width = 10)


