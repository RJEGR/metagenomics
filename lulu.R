#!/usr/bin/env Rscript

# shared=$(ls MAKE_OTUS/*.shared)
# taxonomy=$(ls MAKE_OTUS/*.cons.taxonomy)
# Use: Rscript --vanilla lulu.R downstream/vsearch.match.list $shared $taxonomy
# R version 3.5.0
# Author: Ricardo Gomez-Reyes

# ================
# Defining paths
# ================
path <- paste0(getwd(), "/downstream")
system("mkdir -p downstream")

# ===============
# Check and load package:
# ===============

.git_packages <- c("lulu")
.bioc_packages <- c("biomformat")

.inst <- .git_packages %in% installed.packages()
if(any(!.inst)) {
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools", dep=TRUE, repos='http://cran.us.r-project.org')
    devtools::install_github(.git_packages[!.inst], ask = F)
}


.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.git_packages, .bioc_packages), require, character.only = TRUE)


# # # # # #
# parameters
# # # # # # 

min_r <- 1
min_match <- 98
min_cooccurence <- 0.95 


# # # # # # # #
# Loading data
# # # # # # # #

args = commandArgs(trailingOnly=TRUE)

vsearch.file = args[1]
shared.file = args[2]
taxonomy.file = args[3]

# # # # # # # #
# lulu inputs
# # # # # # # #

#  Match list (vsearch used for all OTUs (Representative seq) with 0.03 % similitud)
matchlist <- read.table(vsearch.file , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

shared <- read.csv(shared.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
otutab <- shared[,-c(1,2,3)]
otutab <- as.data.frame(t(otutab))


# # # # # # # # # # # # # # # #
# define prefix for save resuls:
# # # # # # # # # # # # # # # #
# out.prefix[[1]][-c(1, length(out.prefix[[1]]))]
out.prefix <- strsplit(shared.file, "[.]")
shared.out.prefix <- do.call(paste, 
                             c(as.list(out.prefix[[1]][-c(1, length(out.prefix[[1]]))]), 
                               sep = '.',collapse = ''))

out.prefix <- NULL

out.prefix <- strsplit(taxonomy.file, "[.]")
tax.out.prefix <- do.call(paste, 
                             c(as.list(out.prefix[[1]][-c(1, length(out.prefix[[1]]))]), 
                               sep = '.',collapse = ''))

rds <- paste0(path, "/", "lulu_curated_result.rds")
# options(stringsAsFactors = FALSE)
NO_REUSE = F

if (file.exists(rds) && ! NO_REUSE) {
    print('RESTORING DATA FROM EARLIER ANALYSIS')
    curated_result <- readRDS(rds)
} else {

# # # # # # # # # # #
# Run lulu algorithm
# # # # # # # # # # # 

curated_result <- lulu(otutab, matchlist, minimum_ratio_type = "min", 
                       minimum_ratio = min_r, 
                       minimum_match = min_match, 
                       minimum_relative_cooccurence = min_cooccurence)

# # # # # # # # # # #
# Set results to save
# # # # # # # # # # #

# Save the curated_result object to a file
saveRDS(curated_result, paste0(path, "/", "lulu_curated_result.rds"))
# saveRDS(curated_result, rds)
}

# Restore it later
# curated_result <- readRDS("curated_result.rds")

otu_map <- curated_result$otu_map
parent <- otu_map[otu_map$curated == "parent",]

cat("\n...Number of Processed OTUs are:", length(otu_map$curated), "\n")
cat("Curated in ~real biological OTUs (parents):",table(otu_map$curated)[2], "\n")
cat("And", table(otu_map$curated)[1], "daughters (merged) OTUs\n")

curated_table <- curated_result$curated_table
names(curated_table) <- shared$Group

# # # # # # # # # # # # # #
# recorte de taxonomy file
# # # # # # # # # # # # # #

taxonomy.obj <- read.csv(taxonomy.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
taxonomy <- data.frame(taxonomy.obj[taxonomy.obj$OTU %in% parent$parent_id, ])

tax <- strsplit(taxonomy[,ncol(taxonomy)], ";")
max.rank <- max(lengths(tax))
tax <- sapply(tax, "[", c(1:max.rank)) # Using the max rank assignation to names the taxonomy object
tax <- as.data.frame(t(tax))
tax <- as.data.frame(apply(tax, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)


rank.names <- vector(max.rank, mode="character")

for (i in 1:max.rank) {
  rank.names[i] <- paste("Rank", i, sep="_")
}

colnames(tax) <- rank.names
rownames(tax) <- taxonomy[,1]

# # # # # # # # # # # # 
# Compare files and save
# # # # # # # # # # # #

testing <- identical(rownames(tax), rownames(curated_table))
cat("\n.... Taxonomy and count table has been curated:",testing, "\n")

# # # # # # # #
# Save results
# # # # # # # #

# 1.
write.table(otu_map, 
            file = paste0(path, "/", "lulu", "_",
                          "cooccurence_", min_cooccurence, "_",
                          "min_match_", min_match, "_",
                          "min_ratio_", min_r, "_",
                          "otu_map.csv"),
            sep = " ", quote = FALSE)

# 2.
write.table(curated_table, file = paste0(path,"/" , "lulu.", shared.out.prefix, ".shared"),
            sep = "\t", quote = FALSE)
# 3.
tax.out <- data.frame(OTU = taxonomy[,1],
                      Size = taxonomy[,2], 
                      Taxonomy = apply(tax[,rank.names] , 1 , paste , collapse = ";"))

# 3.1
write.table(tax.out, file = paste0(path, "/", "lulu.", tax.out.prefix, ".taxonomy"),
            sep = "\t", quote = FALSE,
            row.names = FALSE)

# 4
write.table(parent$parent_id, 
            file = paste0(path, "/", "lulu", "_",
                          "cooccurence_", min_cooccurence, "_",
                          "min_match_", min_match, "_",
                          "min_ratio_", min_r, "_",
                          "parent.accnoss"),
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)




# # # # # # # # # 
# Bio integration
# # # # # # # # # 

biom.out <- make_biom(curated_table, 
               observation_metadata = taxonomy,
               sample_metadata = NULL,
               matrix_element_type = "int")

write_biom(biom.out, paste0(path, "/", "lulu", "_",
                            "cooccurence_", min_cooccurence, "_",
                            "min_match_", min_match, "_",
                            "min_ratio_", min_r,
                            ".biom"))


cat("\n Lulu performed ok! ...\n")

quit(save = 'no')