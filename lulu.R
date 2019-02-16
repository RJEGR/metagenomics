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
    devtools::install_github(git_packages[!.inst], ask = F)
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


# # # # # # # # # # # # # # # #
# define prefix for save resuls:
# # # # # # # # # # # # # # # #

out.prefix <- strsplit(shared.file, "[.]")
shared.out.prefix <- do.call(paste, 
                             c(as.list(out.prefix[[1]][-length(out.prefix[[1]])]), 
                               sep = '.',collapse = ''))

out.prefix <- NULL

out.prefix <- strsplit(taxonomy.file, "[.]")
tax.out.prefix <- do.call(paste, 
                             c(as.list(out.prefix[[1]][-length(out.prefix[[1]])]), 
                               sep = '.',collapse = ''))
# # # # # # # #
# Running lulu
# # # # # # # #

#  Match list (vsearch used for all OTUs (Representative seq) with 0.03 % similitud)
matchlist <- read.table(vsearch.file , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

shared <- read.csv(shared.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
otutab <- shared[,-c(1,2,3)]
otutab <- as.data.frame(t(otutab))

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
saveRDS(curated_result, "curated_result.rds")

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
write.table(curated_table, file = paste0(path,"/" ,shared.out.prefix, ".lulu", ".shared"),
            sep = "\t", quote = FALSE)
# 3.
tax.out <- data.frame(OTU = taxonomy[,1],
                      Size = taxonomy[,2], 
                      Taxonomy = tidyr::unite(tax, sep = ";")[,1])

# 3.1
write.table(tax.out, file = paste0(path, "/", tax.out.prefix, ".lulu", ".taxonomy"),
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

quit(save = 'no')

# INCLUIR EN EL RMD ambos archivos (lulu_cons.taxonomy y lulu_curated_table.csv), 

taxonomy2 = observation_metadata(biom.out)
data = as(biom_data(biom.out), "matrix")

# =============

barplot(table(tax[,2]), horiz = TRUE, las=2)

Ph <- as.data.frame(table(tax[,2]))



```{r}
library(lulu)

filtered <- data.frame(otutab[rownames(otutab) %in% parent$parent_id, ])
merged <- data.frame(otutab[!rownames(otutab) %in% parent$parent_id, ])

# compare diversity from both, filtered and raw- count table
# A simple β-diversity measure (average α-diversity divided by γ-diversity) was applied to all uncurated and LULU curated tables (

library(vegan)

# Alfa-div
div.tbl <- data.frame(Id = colnames(filtered),
                      Parent = diversity(t(filtered), index="shannon"),
                      Full = diversity(t(otutab), index="shannon"),
                      Curated = diversity(t(curated_table), index="shannon"))

# Beta-div

div.tbl <- data.frame(Id = colnames(filtered),
                      Parent = betadiver(t(filtered)),
                      Full = betadiver(t(otutab)),
                      Curated = betadiver(t(curated_table)))

Parent.div = betadiver(t(filtered))



library(ggplot2)

# 2 ============good
# Change colors
ggplot(div.tbl, aes(x=Parent, y=Curated, color=Curated)) +
  geom_point(aes(size = Curated)) + geom_rug()
# rownames(div.tbl) <- sample_data(phyloseq)$`Estación`

# 3 ============

p <- ggplot(div.tbl, aes(Parent, Curated, shape = factor(cyl)))
p + geom_point(aes(colour = factor(cyl)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)

div.plot <- reshape2::melt(div.tbl)
# 4 ======
p <- ggplot(div.plot, aes(x=variable, 
                          y=value, 
                          color=variable)) + 
  geom_boxplot(notch=TRUE) + coord_flip()

# Use grey scale
p + scale_color_grey() + theme_classic()
# p + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))





# include and subset taxonomy file, then save the file.taxonmy newer and output to the dowstream folder

# 

ggplot(otu_map, aes(x=spread, y=log2(total))) + 
  geom_point(aes(color=curated, shape=curated)) +
  scale_color_manual(values=c("#E69F00", "#999999")) +
  facet_wrap( ~ curated)
# or

# density plot!
#Figure 7. Distribution of best matches for OTUs on GenBank.
# 
#Density distribution of the best reference database match for all OTUs (percent
# identity (%) of best matching reference sequence on GenBank) is plotted as a
# violin plot. Red bars represent OTUs discarded by the LULU algorithm, and blue
# bars represent OTUs retained/curated by the algorithm

#  ======= Redundancy plots
# taxonomic redundancy (percentage of OTUs with redundant taxonomic annotation)
# re-incorporar taxones
taxa1 <- as.matrix(tax_table(phyloseq))
taxa <- data.frame(taxa1[rownames(taxa1) %in% parent$parent_id, ]) # time demand!
Rank <- as.vector(taxa1[,1])
tbl.plot <- cbind(otu_map, Rank)

ggplot(tbl.plot, aes(x=spread, y=log2(total))) + 
  geom_point(aes(color=Rank, shape=curated)) +
  #scale_color_manual(values=c("#E69F00", "#999999")) +
  facet_wrap( ~ Rank)
# based on the plot above
tail(tbl.plot[tbl.plot$curated == 'merged',])
tbl.plot[tbl.plot$Rank == 'Plantae',]
table(tbl.plot[tbl.plot$Rank == 'Plantae',][,'curated'])

links <- data.frame(node1 = tbl.plot[tbl.plot$Rank == 'Plantae',]['parent_id'],
                    node2 = rownames(tbl.plot[tbl.plot$Rank == 'Plantae',]))

# and net

library(igraph)
vertex <- tbl.plot[tbl.plot$Rank == 'Plantae',]
vertex <- cbind(node1 = rownames(vertex), vertex)
net <- graph_from_data_frame(d = links, # nodes interaction
                             vertices = vertex, # attribudes per node
                             directed=TRUE)

plot(degree_distribution(net),
     #log="xy",
     col="sienna4",
     xlab = "Number of connections",
     ylab = "Density")


net <- igraph::simplify(net, remove.multiple = F, remove.loops = T) 

plot(net, layout=layout_with_kk,
     edge.arrow.size=.2, # reduce arrow size
     vertex.size=5,
     vertex.label=NA,
     edge.curved=.1)
```
