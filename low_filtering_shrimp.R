rm(list = ls())

# lets EDA 

library(tidyverse)
library(ANCOMBC)
library(phyloseq)
library(microbiome)

dir <- '~/Documents/Shrimp_Estefany/'

ps <- readRDS(paste0(dir, "ps_normalized.rds"))

# test outliers

# Look for low perfroming samples

ps %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  tax_glom(., "Family", NArm = FALSE) -> ps2

qplot(colSums(otu_table(ps)),bins=30) + 
  xlab("Logged counts-per-sample")
