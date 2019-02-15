#!/usr/bin/env Rscript

setwd(dir="/Users/cigom/metagenomics/LULU")
# for i in $(seq 0.010 0.002 0.052); do Rscript --vanilla tax.summary.R $i; done
library(lulu)
# :::::::::::::::::: loading args

args = commandArgs()
cutoff = args[7]
path <- paste0("otutab/", "cuttoff_", cutoff,"/")
path2 <- paste0("taxonomy/", "otus_", cutoff, "/")
file1 <- dir(path = path, pattern = ".shared")
file2 <- dir( path = path2, pattern = ".cons.taxonomy")

file <- paste0(path, file1[1])
print(file)

shared <- read.csv(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
otutab <- shared[,-c(1,2,3)]
otutab <- as.data.frame(t(otutab))

names <- strsplit(shared[,2], "_")
names <- sapply(names, "[", c(3))

colnames(otutab) <- names

otutab$cutoff <- cutoff
# ::: parse classification also
file <- NULL
file <- paste0(path2, file2[1])
print(file)

constaxa <- read.csv(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
# constaxa <- read.csv("MAKE_OTUS/cigom.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.cons.taxonomy", sep="\t", header=TRUE, stringsAsFactors=FALSE)

tax <- strsplit(constaxa[,3], ";")
tax <- sapply(tax, "[", c(1:7))
tax <- as.data.frame(t(tax))

Ranks <- c("Reino", "Filo", "Clase", "Orden", 
                   			 "Familia", "Genero", "Especie")

T <- as.data.frame(apply(tax, 2, 
    function(x) gsub("\\(.*$", "",  x, perl=TRUE)),
     stringsAsFactors = F)

colnames(T) <- Ranks                              
# ::::::::::::;

save <- cbind(otutab, T)

write.table(save, file = paste0(cutoff,".shared.taxonomy", ".csv"),
                sep = " ", quote = FALSE, row.names = TRUE, col.names = TRUE)


quit(save = "no")


