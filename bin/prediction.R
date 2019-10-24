#!/usr/bin/env Rscript

library("tidyverse")

# parameters -------------------------------------------------------------------
# 1: ps_copyadj
# 2: ps_copyadj_val

args <- commandArgs(trailingOnly = TRUE)

ps_data <- readRDS(args[[1]])
ps_data_val <- readRDS(args[[2]])

phyloseq::sample_data(ps_data)$batch <- "first_run"
phyloseq::sample_data(ps_data_val)$batch <- "second_run"

# don't want phy tree
ps_no_tree <- phyloseq::phyloseq(
  phyloseq::sample_data(ps_data),
  phyloseq::otu_table(ps_data),
  phyloseq::tax_table(ps_data)
)
ps_val_no_tree <-
  phyloseq::phyloseq(
    phyloseq::sample_data(ps_data_val),
    phyloseq::otu_table(ps_data_val),
    phyloseq::tax_table(ps_data_val)
  )
ps_all <- phyloseq::merge_phyloseq(ps_no_tree, ps_val_no_tree)
ps_all <- phyloseq::subset_samples(ps_all, smoking != "Uncertain")

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

prevalence_threshold <- 0.05 * phyloseq::nsamples(ps_all)
keep_taxa <- names(prevdf)[(prevdf >= prevalence_threshold)]
dep_data <- phyloseq::prune_taxa(keep_taxa, ps_all)

# convert counts to relative abundance
rel_dep_data <-
  phyloseq::transform_sample_counts(dep_data, function(x)
    x / sum(x))

dep_tab <-
  data.matrix(unclass(t(phyloseq::otu_table(rel_dep_data))))
dep_tab <-
  tibble::as_tibble(dep_tab[,-caret::nearZeroVar(dep_tab)])
dep_samp_df <-
  data.frame(
    class = phyloseq::sample_data(rel_dep_data)$cohort,
    smoking = phyloseq::sample_data(rel_dep_data)$smoking,
    batch = phyloseq::sample_data(rel_dep_data)$batch
  )
dep_df <- bind_cols(dep_samp_df, dep_tab) %>%
  mutate(smoking = factor(smoking, levels = c("Never", "Daily"))) %>%
  filter(!is.na(smoking))

dep_df %>%
  filter(batch == "first_run") %>%
  select(-class, -smoking, -batch) -> train
dep_df %>%
  filter(batch == "first_run") %>%
  select(class, smoking) %>%
  mutate(smoking = as.numeric(smoking))-> train_meta


dep_df %>%
  filter(batch == "second_run") %>%
  select(-class,-smoking,-batch) -> test
dep_df %>%
  filter(batch == "second_run") %>%
  select(class, smoking) %>%
  mutate(smoking = as.numeric(smoking)) -> test_meta

train <- scale(train)
test <- scale(
  test,
  center = attr(train, "scaled:center"),
  scale = attr(train, "scaled:scale"))

# coerce to data frame after scaling
# (don't do this before or scaling attributes get lost)
train_meta %>%
  bind_cols(data.frame(train)) %>%
  select(class, everything()) -> train_df

test_meta %>%
  bind_cols(data.frame(test)) %>%
  select(class, everything()) -> test_df

set.seed(0451)
model <-
  OmicsMarkeR::fs.ensembl.stability(
    as.matrix(train_df[,-1]),
    # drop class
    train_df[, 1],
    # class labels
    method = c("svm", "rf"),
    k = 5,
    # number of bootstrapped iterations,
    p = 0.7,
    f = 10,
    metric = "Kappa",
    allowParallel = FALSE # parallel messes up seed
  )

sink("prediction.txt")
OmicsMarkeR::performance.metrics(model)

preds <- OmicsMarkeR::predictNewClasses(model, method = "rf", train_df, test_df[, -1])
caret::confusionMatrix(preds$predictedClass, test_df[, 1], positive = "Depression")
sink()