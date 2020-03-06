#!/usr/bin/env Rscript

library("tidyverse")

# parameters -------------------------------------------------------------------
# 1: lefse results file

args <- commandArgs(trailingOnly = TRUE)

swr = function(string, nwrap=40) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

sws <- tibble::as_tibble(read.csv(args[[1]]))
sws <- sws %>% group_by(cohort) %>%
  mutate(lda = ifelse(cohort == "Depression", 0-lda, lda)) %>%
  mutate(description = reorder(desc, lda)) %>%
  mutate(superclass = swr(superclass)) %>%
  mutate(description = swr(description))

ggplot(sws, aes(x = description, y = lda, label = description, fill = cohort)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + 
  coord_flip() +
  xlab("") + 
  ylab("LDA SCORE (log 10)") + 
  scale_color_brewer(palette="Set1", direction = 1) +
  facet_wrap(~superclass, scales="free")

# levels(factor(sws$superclass))
# nah pick top 4 superclasses
keep_lvls <- c("Cofactor, Prosthetic Group, Electron\nCarrier, and Vitamin Biosynthesis",
          "Carbohydrate Biosynthesis",
          "Nucleoside and Nucleotide Biosynthesis",
          "Carbohydrate Biosynthesis")
sws %>%
  mutate(superclass = forcats::fct_other(superclass, keep_lvls, other_level = "Other")
) -> sws_other

ggplot(sws_other, aes(x = description, y = lda, label = description, fill = cohort, width = 0.75)) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  coord_flip() +
  xlab("") + 
  ylab("LDA SCORE (log 10)") + 
  scale_color_brewer(palette="Set1", direction = 1) +
  facet_wrap(~superclass, scales="free") + 
  theme(legend.position = "none") +
  theme(text = element_text(size=22))

ggsave("picrust.pdf", device = "pdf", width = 25, height = 40)
ggsave("picrust.svg", device = "svg", width = 25, height = 40)

# save each facet individually
# Create a separate plot for each value of cyl, and store each plot in a list
p.list = lapply(sort(unique(sws_other$superclass)), function(i) {
  ggplot(sws_other[sws_other$superclass==i,], aes(x = description, y = lda, label = description, fill = cohort, width = 0.5)) +
    geom_bar(stat = "identity") + 
    theme_bw() + 
    coord_flip() +
    xlab("") + 
    ylab("LDA SCORE (log 10)") + 
    scale_color_brewer(palette="Set1", direction = 1) +
    facet_wrap(~superclass, scales="free") + 
    theme(legend.position = "none") +
    theme(text = element_text(size=14)) 
})

# save manually in RStudio!
