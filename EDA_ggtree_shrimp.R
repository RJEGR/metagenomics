# wrangling pos-clustering results
# https://yulab-smu.top/treedata-book/chapter6.html
.cran_packages <- c("dplyr")
.bioc_packages <- c("Biostrings", "DECIPHER", "phangorn")

sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

library(Biostrings)

fasta_file <- paste0(dir, "dna-sequences.fasta")
seqs <- Biostrings::readDNAStringSet(fasta_file)

# Construct phylogenetic tree

alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose = FALSE)

# The phangorn R package is then used to construct a phylogenetic tree. 
# follow lines construct a neighbor-joining tree, 
# and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree 
# using the neighbor-joining tree as a starting point

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)

treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

save(fitGTR,fit, treeNJ, dm,phangAlign, file = paste0(dir, "dna-sequences-fitGTR.RData"))

detach("package:phangorn", unload=TRUE)

library(tidyverse)

tree <- fitGTR$tree

load(paste0(dir, "lulu_curated_result.Rdata"))
load(paste0(dir, "Estefy_rarefied_table_w_taxonomy.RData"))

otu_map <- curated_result$otu_map %>% as_tibble(rownames = "label")

names(obj)[1] <- 'label'

obj %>% select_if(is.character) -> tax
obj %>% select_if(negate(is.character))-> ab

library(ggtree)

# BiocManager::install("ggtree")


x <- as_data_frame(tree)
y <- full_join(x, otu_map, by = 'label') %>% full_join(tax) %>% arrange(curated)

y %>% pull(curated) %>% unique()

rn <- y$parent_id

group <- list(merged = rn[which(y$curated %in% "merged")], 
     parent = rn[which(y$curated %in% "parent")])

# el objeto y es el bueno para explorarlo, como? interesante!

tre <- NULL
tre <- y %>% tidytree::as.treedata()

# library(RColorBrewer)
# 
# colourCount = length(unique(otu_map$curated))
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggtree(tre) + layout_inward_circular()

# 

ggtree(tre, layout = 'circular', branch.length='none') %>%
  groupOTU(., group, 'Curated') + aes(color=curated) +
  theme(legend.position="right") +
  scale_color_manual(values = c("red", "blue"),
                     breaks = c("merged", "parent"))
