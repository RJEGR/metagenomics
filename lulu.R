#!/usr/bin/env Rscript

# module load R-3.5.0
# shared=$(ls MAKE_OTUS/*.shared)
# taxonomy=$(ls MAKE_OTUS/*.cons.taxonomy)
# Use: Rscript --vanilla lulu.R downstream/vsearch.lulu_in.list.txt $shared $taxonomy
# R version 3.5.0
# Author: Ricardo Gomez-Reyes

# ================
# Defining paths
# ================
path <- getwd()
# path <- '~/metagenomics/run13_18S/downstream/'

# ===============
# Loading package:
# ===============

if (!require('devtools')) {
  install.packages("devtools", dep=TRUE, repos='http://cran.us.r-project.org')
} else   
  if (!require('lulu')) {
    devtools::install_github("tobiasgf/lulu")
  }

# # # # # # # #
# Loading data
# # # # # # # #


args = commandArgs(trailingOnly=TRUE)

vsearch.file = args[1]
shared.file = args[2]
taxonomy.file = args[3]

if (is.na(args[4])) {
  Rank <- c("Reino", "Filo", "Clase", "Orden", "Familia", "Especie")
} else
  Rank <- c(args[4])


# # # # # # # #
# Running lulu
# # # # # # # #

#  Match list (used vsearch for all OTUs (Representative seq) with 0.03 % similitud)
matchlist <- read.table(vsearch.file , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

shared <- read.csv(shared.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
otutab <- shared[,-c(1,2,3)]
otutab <- as.data.frame(t(otutab))

# # # # # # # # # # #
# Run lulu algorithm
# # # # # # # # # # # 

curated_result <- lulu(otutab, matchlist, minimum_ratio_type = "min", 
                       minimum_ratio = 1, 
                       minimum_match = 98, 
                       minimum_relative_cooccurence = 0.95)


# 
# otu_map <- read.csv(paste0(path, "/lulu_otu_map", ".csv"), sep = ' ')
# curated_table <- read.csv(paste0(path,"/lulu_curated_table", ".csv"), sep = ' ')

library(dplyr)

parent <- filter(otu_map, curated == "parent")
filtered <- data.frame(otutab[rownames(otutab) %in% parent$parent_id, ])
merged <- data.frame(otutab[!rownames(otutab) %in% parent$parent_id, ])

# # # # # # # # # # # # # #
# recorte de taxonomy file
# # # # # # # # # # # # # #

# taxonomy.file = 'cigom.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.cons.taxonomy'

taxonomy.obj <- read.csv(taxonomy.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
taxonomy <- data.frame(taxonomy.obj[taxonomy.obj$OTU %in% parent$parent_id, ])

tax <- strsplit(taxonomy[,ncol(taxonomy)], ";")
tax <- sapply(tax, "[", c(1:length(Rank))) # solo se toman los niveles de Rank pero se puede tomar el resto
tax <- as.data.frame(t(tax))
tax <- as.data.frame(apply(tax, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)


colnames(tax) <- Rank
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

otu_map <- curated_result$otu_map
write.table(otu_map, 
            file = paste0(path, "/downstream/lulu_otu_map", ".csv"),
            sep = " ", quote = FALSE)

# 2.

curated_table <- curated_result$curated_table
write.table(curated_table, file = paste0(path,"/downstream/lulu_list.shared", ".csv"),
            sep = " ", quote = FALSE)
# 3.

tax.out <- data.frame(OTU = taxonomy[,1],
                      Size = taxonomy[,2], 
                      Taxonomy = tidyr::unite(tax, sep = ";")[,1])

write.table(tax.out, file = paste0(path,"/downstream/lulu_cons.taxonomy"),
            sep = "\t", quote = FALSE,
            row.names = FALSE)

# # # # # # # # # 
# Bio integration
# # # # # # # # # 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomformat", version = "3.8")

library(biomformat)

make_biom()

#MAKE_OTUS/*.cons.taxonomy # este es el original
#downstream/*.cons.taxonomy # el archivo es el mismo , pero contiene el Size del numero de la secuencia representativa





quit(save = 'no')

# =============

barplot(table(tax[,2]), horiz = TRUE, las=2)

Ph <- as.data.frame(table(tax[,2]))



```{r LULU}
library(lulu)

# Loading match list (used vsearch for all OTUs (Representative seq) with 0.03 % similitud)

# downstream <- paste0(path, "/downstream")
vsearch.file <- paste0(path,"/downstream", "/vsearch.lulu_in.list.txt")
matchlist <- read.table(vsearch.file , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

otutab <- as.data.frame(otu_table(phyloseq), header=TRUE, stringsAsFactors=FALSE)
# otutab <- as.data.frame(t(otutab))

#names <- strsplit(names(shared), "_")
#names <- sapply(names, "[", c(3))
#names(shared) <- names

curated_result <- lulu(otutab, matchlist, minimum_ratio_type = "min", 
                       minimum_ratio = 1, 
                       minimum_match = 98, 
                       minimum_relative_cooccurence = 0.95)

# And update phylose object

otu_map <- curated_result$otu_map

write.table(otu_map, 
            file = paste0(path, "/lulu_otu_map", ".csv"),
            sep = " ", quote = FALSE)


table(otu_map$curated)
# merged parent (run13_18S)
#  4221  22101


# 
# :::::: also save the newer table

curated_table <- curated_result$curated_table

write.table(curated_table, file = paste0(path,"/downstream/lulu_curated_table", ".csv"),
            sep = " ", quote = FALSE)

# 

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
