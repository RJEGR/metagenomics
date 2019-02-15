#!/usr/bin/env Rscript

setwd(dir="/Users/cigom/metagenomics/LULU")
# for i in $(seq 0.010 0.002 0.052); do Rscript --vanilla lulu_cuttoff.R $i; done
library(lulu)
# :::::::::::::::::: loading args

args = commandArgs()
cutoff = args[7]
# cutoff = 0.03
path <- paste0("otutab/", "cuttoff_",cutoff,"/")
path2 <- paste0("vsearch/", "vsearch.")
file1 <- dir(path = path, pattern = ".shared")

# print(paste0(path, file1, " and ", path2, cutoff, ".list.txt"))

# ::::::::::::::::::::::::: Loading the per-sample OTU matrix
file <- paste0(path, file1[1])
print(file)

shared <- read.csv(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# shared <- read.csv("MAKE_OTUS/cigom.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.shared", sep="\t", header=TRUE, stringsAsFactors=FALSE)
otutab <- shared[,-c(1,2,3)]
otutab <- as.data.frame(t(otutab))

names <- strsplit(shared[,2], "_")
names <- sapply(names, "[", c(3))

colnames(otutab) <- names


# :::::::::::::::::::
# ::::::::::::::::::: LOADING MATCH LIST
# :::::::::::::::::::
file2 <- paste0(path2, cutoff, ".list.txt")
matchlist <- read.table(file2 , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)
# matchlist <- read.table("vsearch_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

#:::::::::::::::::::::::: Run lulu

curated_result <- lulu(otutab, matchlist, minimum_ratio_type = "min", 
                            minimum_ratio = 1, 
                            minimum_match = 98, 
                            minimum_relative_cooccurence = 0.95)

otu_map <- curated_result$otu_map
otu_map$cutoff <- cutoff

write.table(otu_map, file = paste0(cutoff,"_otu_map", ".csv"),
                sep = " ", quote = FALSE)

# :::::: also save the newer table

curated_table <- curated_result$curated_table
curated_table$cutoff <- cutoff

write.table(curated_table, file = paste0(cutoff,"_curated_table", ".csv"),
                sep = " ", quote = FALSE)

quit(save = 'no')

