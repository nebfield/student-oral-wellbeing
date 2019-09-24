#!/usr/bin/env Rscript
# TODO!
library("tidyverse")

ps_data <- readRDS("~/Downloads/ps_copyadj.rds")
second_batch <- 
ps_data <- subset_samples(ps_data, smoking != "Uncertain")

# step one: filter low-abundance bugs
prevdf <- apply(
  X = phyloseq::otu_table(ps_data),
  MARGIN = ifelse(
    phyloseq::taxa_are_rows(ps_data),
    yes = 1,
    no = 2
  ),
  FUN = function(x) {
    sum(x > 0)
  }
)

prevalence_threshold <- 0.05 * phyloseq::nsamples(ps_data)
keep_taxa <- names(prevdf)[(prevdf >= prevalence_threshold)]
dep_data <- phyloseq::prune_taxa(keep_taxa, ps_data) 

# convert counts to relative abundance
rel_dep_data <- phyloseq::transform_sample_counts(dep_data, function(x) x/sum(x))

dep_df <-
  data.frame(
    class = phyloseq::sample_data(rel_dep_data)$cohort,
    smoking = phyloseq::sample_data(rel_dep_data)$smoking,
    data.matrix(unclass(t(phyloseq::otu_table(rel_dep_data))))
  )

scale_data <- function(data_bal) {
  data_bal <- data_bal[, -caret::nearZeroVar(data_bal)] 
  
  x <- data_bal %>% 
    dplyr::select(-class,-smoking) # train + test
  
  x <- scale(x)
  
  # coerce to data frame after scaling
  # (don't do this before or scaling attributes get lost)
  x <- data.frame(x) 
  
  # add smoking data back
  x$smoking <- as.numeric(data_bal$smoking)
  
  # add class labels back
  x$class <- data_bal$class
  
  # make sure class label is first column 
  return(x %>% select(class, everything()))
}

dep <- scale_data(data = dep_df)

set.seed(0451)
model <-
  OmicsMarkeR::fs.ensembl.stability(
    as.matrix(dep[, -1]),
    # drop class
    dep[, 1],
    # class labels
    method = c("svm", "rf"),
    k = 5, # number of bootstrapped iterations,
    p = 0.7,
    f = 10,
    metric = "Kappa",
    allowParallel = FALSE # parallel messes up seed
  )

OmicsMarkeR::performance.metrics(model)
