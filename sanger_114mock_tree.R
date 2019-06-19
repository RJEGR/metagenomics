#!/usr/bin/env Rscript

# ===+===
# Treeing
# ===+====
cat("\n**********\n**********\nConstruct a neighbor-joining and then fit a GTR+G+I maximum likelihood tree\n**********\n**********\n\n\n")

# =====================
# Loading data results
# =====================
rm(list = ls())

args = commandArgs(trailingOnly=TRUE)

path <- '/Users/cigom/metagenomics/COI/run012/mock_P68/mock_parameter_definition'
# path_files <- '/Users/cigom/metagenomics/sanger_assign_db/'
# fasta.file <- paste0(path_files, 'OMEGA_A_1e_120_ASVs.fasta')
# out_prefix <- 'OMEGA_A_1e_120_ASVs'

# bold_sp <- 'OMEGA_A_1e_120_ASVs.BOLD_public_species.wang.taxonomy'
# taxonomy.file <- paste0(path_files, bold_sp)
setwd(path)
# file.name <- args[1]
fasta.file <- 'ictio_coi_sanger114.fasta'
out_prefix <- 'ictio_coi_sanger114'

taxonomy.file <- 'ictio_coi_sanger114.tax'

if (is.na(args[2])) {
  threads <- NULL
} else
  threads = as.numeric(args[2])

# ================
## Checking and Load packages ----
# ================

.cran_packages <- c("dplyr")
.bioc_packages <- c("Biostrings", "DECIPHER", "phangorn", "ggtree", "treeio")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F, version = "3.8")
}


# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# ================
# Outputs in `pwd`:
# ================
dir = getwd()
# dir = path_files
out_path <- file.path(dir, "")
system(command = paste0("mkdir -p ", out_path), intern = F)

# ===========================

# Construct phylogenetic tree
seqs <- readDNAStringSet(fasta.file)

# seqs <- seqs[names(seqs) %in% rownames(boots),]

seqs.width <- width(seqs)
names_seqs <- names(seqs)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, processors = threads)

# The phangorn R package is then used to construct a phylogenetic tree. 
# follow lines construct a neighbor-joining tree, 
# and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree 
# using the neighbor-joining tree as a starting point

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)

# Distance based methods are very fast and we will use the UPGMA and NJ tree as starting trees for the maximum parsimony and maximum likelihood analyses.

treeNJ <- NJ(dm) # Note, tip order != sequence order

# treeUPGMA <- upgma(dm)
# layout(matrix(c(1,2), 2, 1), height=c(1,2))
# par(mar = c(0,0,2,0)+ 0.1)
# plot(treeUPGMA, main="UPGMA", show.tip.label = FALSE)
# plot(treeNJ, "unrooted", main="NJ", show.tip.label = FALSE)
# Parsinomy score means the number of changes which are at least necessary to describe the data for a given tree
# parsimony(treeUPGMA, phangAlign) # 5136
# parsimony(treeNJ, phangAlign) # 5037


fit = pml(treeNJ, data=phangAlign)

fitGTR <- update(fit, k=4, inv=0.2)

fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)

# save.image(file = paste0(getwd(),"/", out_prefix,".results.RData"))

tree <- fitGTR$tree

# tree$tip.seq <- tree$tip.label
# tree$tip.label <- names_seqs

saveRDS(tree, paste0(getwd(),"/", out_prefix, "_fitGTR-tree.rds"))
# tree <- readRDS(paste0(getwd(),"/", out_prefix, "_fitGTR-tree.rds"))

cat("\nOutputs saved at\n", paste0(getwd(),"/", out_prefix, "_fitGTR-tree.rds"),"\n")



# quit(save='no')

cat("\n 1. Processing taxonomy... \n")

taxonomy.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) # Using the max rank assignation to names the taxonomy object
taxonomy <- as.data.frame(t(taxonomy))
rownames(taxonomy) <- taxonomy.obj[,1]
tax <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)
tax <- mutate_all(data.frame(tax), funs(str_replace_all(., c("_unclassified"="", "Unclassified"="", "_"=" "))))


rank.names <- vector(max.rank, mode="character")

for (i in 1:max.rank) {
  rank.names[i] <- paste("Rank", i, sep="_")
}

rank.names <- c("Reino", "Filo", "Clase", "Orden", "Familia", "Genero", "Especie")
# rank.names <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax) <- rank.names
rownames(tax) <- taxonomy.obj$V1

tax <- as_tibble(data.frame(label = rownames(tax), tax, row.names =  NULL))

data.frame(taxa=apply(tax, 2, function(x) length(table(x)) ))

# Rank_1 Rank_2 Rank_3 Rank_4 Rank_5 Rank_6 Rank_7 
#   1      1      1     10     34     57     69 

# GGTREE

tbl <- as_tibble(read.csv("/Users/cigom/metagenomics/COI/X04_SangerIctio/ListaPool114_sanger.csv", sep=",", stringsAsFactors = FALSE, header=TRUE))

names(tbl)[1] <- "label"
names(tbl)[3] <- "Species"

library(ggplot2)
library(RColorBrewer)

tree$tip.label <- names_seqs

x <- as_data_frame(tree)

x$label <- names_seqs

x %>% full_join(tbl, by = 'label') %>%
      full_join(tax, by = 'label') -> y

x %>% 
  full_join(tax, by = 'label') %>%
  anti_join(tbl, by = 'Species') -> y

tree <- NULL
tree <- y %>% tidytree::as.treedata()

# y %>% filter(!is.na(y[,rank.names[7]])) 
colourCount = nrow(unique(y[,'Order']))
getPalette = colorRampPalette(brewer.pal(colourCount, "Paired"))

plot <- ggtree(tree, branch.length='none', layout='slanted', aes(color=Order,  na.rm = TRUE), size=1) +
  theme(legend.position="top") +
  scale_color_manual(values = c(getPalette(colourCount)), na.value = "grey", guide = guide_legend(ncol=round(colourCount / 3))) +
  #geom_tiplab(size=4, aes(label=paste0('italic(', label, ')~bolditalic(', Query.ID, ')~', Best.ID)), parse=FALSE)
  xlim(0, 40) +
  geom_tiplab(size=2.5, aes(label=Best.ID), offset=0.2) + # color = 'black'
  geom_nodepoint(alpha=1/4, size=2.5) +
  labs(caption = "GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree")


png(paste0(getwd(),"/", out_prefix, "_fitGTR-tree-slated_plot.png"), units="px", width=2500, height=4000, res=400)
plot(plot, col=adjustcolor("black", alpha=0.2))
dev.off()

# 2circular
plot <- NULL

plot <- ggtree(tree, layout='circular', aes(color=Rank_4), na.rm = TRUE, size=1) +
  theme(legend.position="top") +
  scale_color_manual(values = c(getPalette(colourCount)), 
                     na.value = "grey", guide = guide_legend(ncol=round(colourCount / 2))) +
  geom_tiplab2(size=2.5, aes(label=Best.ID), offset=0.2) +
  labs(caption = "GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree")
  
png(paste0(getwd(),"/", out_prefix, "_fitGTR-tree-circular_plot.png"), units="px", width=4000, height=4000, res=400)
plot(plot, col=adjustcolor("black", alpha=0.2))
dev.off()

# rectangular
plot <- NULL
plot <- ggtree(tree, branch.length='none', layout='rectangular', aes(color=Rank_4,  na.rm = TRUE), size=1) +
  theme(legend.position="top") +
  scale_color_manual(values = c(getPalette(colourCount)), na.value = "grey", guide = guide_legend(ncol=round(colourCount / 3))) +
  xlim(0, 40) +
  geom_tiplab(size=2.5, aes(label=Best.ID), offset=0.2) + # color = 'black'
  geom_nodepoint(alpha=1/4, size=2.5) +
  labs(caption = "GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree")

png(paste0(getwd(),"/", out_prefix, "_fitGTR-tree-rect_plot.png"), units="px", width=2500, height=4200, res=400)
plot(plot, col=adjustcolor("black", alpha=0.2))
dev.off()

# plot nodes
# y %>% filter(is.na(y[,7])) 
ggtree(tree) +
  geom_hilight(node=c(116, 122), fill="darkgreen", alpha=.6) +
  geom_hilight(node=c(137,118), fill="steelblue", alpha=.6)

           






