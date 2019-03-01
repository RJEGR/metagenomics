#!/usr/bin/env Rscript

# Use: Rscript --vanilla luluVis.R downstream/vsearch.match.list $shared $taxonomy
# R version 3.5.0
# Author: Ricardo Gomez-Reyes

# ===============
# Check and load package:
# ===============

.bioc_packages <- c("biomformat")


.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

sapply(c(.bioc_packages), require, character.only = TRUE)


# setdir
setwd('/Users/cigom/metagenomics/run13_18S/downstream')

# 1.
curated_result <- readRDS("curated_result.rds")
otu_map <- curated_result$otu_map


# 2. 
taxonomy.file = 'cigom.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.cons.taxonomy'
taxonomy.obj <- read.csv(taxonomy.file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) # Using the max rank assignation to names the taxonomy object
taxonomy <- as.data.frame(t(taxonomy))
rownames(taxonomy) <- taxonomy.obj[,1]

tax <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)

for(i in ncol(tax):1) {
    # it will replace 'unknown_unclassified' terms for 'Unknown' term
    tax[,i] <- sub("unknown", "Unknown", tax[,i])
    tax[,i] <- sub("unknown.*", "Unknown", tax[,i])
#   Sub taxa with *nclassified with Undetermined term
    tax[,i] <- sub(".*nclassified.*", 
    "Undetermined", tax[,i], perl=TRUE)
#   Instead of Undetermined classified , use NA
    #tax[,i][grep("Undetermined|Unknown|NA",
    tax[,i][grep("Undetermined|NA",
                        tax[,i], perl=TRUE)] <- NA
}

# ===

# variales
rank <- 1

rank.names <- vector(max.rank, mode="character")

for (i in 1:max.rank) {
  rank.names[i] <- paste("Rank", i, sep="_")
}

colnames(tax) <- rank.names
rownames(tax) <- taxonomy.obj[,1]

# # # # # # #
# Visualizing
# # # # # # #
library(ggplot2)


# 1. Distribucion de las asignaciones (ie. asigaciones y NAs) que 
#hay entre los OTUs parentales y los que se absorben

lineage <- data.frame(lineage = tax[,rank], curated = otu_map$curated)
lineage.stats <- aggregate(lineage[,2], by=list(lineage[,1]), FUN = table)
datavis <- melt(data.frame(lineage = lineage.stats[,1], lineage.stats[,2]))


datavis$Pct <- round(datavis$value / sum(datavis$value) * 100)

ggplot(data=datavis, aes(x = lineage, y = Pct, fill = variable))+
        geom_bar(stat="identity", colour="black",size=0.05, position = position_stack(reverse = TRUE))+
        #coord_polar()+
        coord_flip()+
        #geom_text(aes(label=value), size=3, vjust=-0.2) +
        geom_text(data = subset(datavis, value > 0), aes(label=value), size=3, vjust=-0.2) +
        scale_fill_brewer(palette="Greens") +
        theme_classic() + ylim(0,100) +
        labs(
            y = "% de amplicones",
            x = paste0("Nivel taxonómico ", " (",names(parent.tax)[rank],")"), 
            title="Redundancia en la asignación taxonómica", 
            subtitle=paste0("OTUs absorbidos.\nTotal amplicones = ", nrow(tax), 
                                    " (Parent = ", nrow(parent.tax), ")",
                                    " (Merged = ", nrow(merged.tax), ")"))

# done!

# 2. Como varia la diversidad (alfa, beta (T.G Frosvel et al y dbotu3 implementan B-diversidad) 
# a lo largo de las muestras (boxplot, se refleja la inflacion de la estimacion de la diversidad)

library(vegan)

dim(original_table <- curated_result$original_table)
dim(curated_table <- curated_result$curated_table)
# Alfa-div

library(reshape2)

alfa.tbl <- melt(data.frame(
                      original_table = diversity(t(original_table), index="shannon"),
                      curated_table = diversity(t(curated_table), index="shannon")))

# datavis <- reshape2::melt(alfa.tbl)

# compare diversity from both, filtered and raw- count table
# A simple β-diversity measure (average α-diversity divided by γ-diversity) was applied to all uncurated and LULU curated tables (

richness <- colSums(curated_table)
obs_beta <- nrow(curated_table)/mean(richness)

ggplot(alfa.tbl, aes(x=variable, 
                          y=value, 
                          color=variable)) + 
  geom_boxplot(notch=TRUE) + coord_flip() +
  scale_color_grey() + theme_classic() +
  labs(x = "",
       y = "Alfa diversity")

# grey = '#999999'; yellow = '#E69F00'; blue = '#56B4E9'
# p + scale_color_manual(values=c(grey, yellow))

original_table_ = betadiver(t(original_table))
curated_table_ = betadiver(t(curated_table))
plot(original_table_)

# 3. Distribucion de OTUs parentales/absorbidos vs su dispersion a lo largo de las muestras (scatterplot colored)

# include and subset taxonomy file, then save the file.taxonmy newer and output to the dowstream folder

# este grafico no dice nada, 
ggplot(otu_map, aes(x=spread, y=log2(total))) + 
  geom_point(aes(color=curated, shape=curated), alpha = 0.7) +
  scale_color_manual(values=c("#E69F00", "#999999")) +
  xlim(0,max(otu_map$spread)) +
  theme_classic() +
  facet_wrap( ~ curated) 

# Amplicon spread through the samples
#xlim <- c(max(otu_map$spread))

ggplot(otu_map, aes(x=spread, color = curated, fill = curated)) +
  geom_histogram(aes(y=..count..), bins = 100, position = "dodge") +
  scale_color_grey() + scale_fill_grey() +
  theme_classic() +
  labs(
            x = "n muestras", 
            y = "Amplicones",
            title="Dispersión de amplicones a lo largo de las muestras", 
            subtitle=paste0( "Total amplicones = ", nrow(tax), 
                                    " (Parent = ", nrow(parent.tax), ";",
                                    " Merged = ", nrow(merged.tax), ")", 
                                    " a lo largo de ", max(otu_map$spread), " muestras"))

# 4. Red de OTUs absorbidos por parentales coloreado su asignacion taxonomica. 
# (quiza solo se presente el numero de asignaciones NA que fueron abosrbidas en un parental)

links <- data.frame(node1 = rownames(otu_map),
                    node2 = otu_map$parent_id,
                    index = otu_map$curated,
                    s_samples = otu_map$spread,
                    lineage = tax[,rank])

library(tidygraph)
library(ggraph)

# install.packages('networkD3')
routes_tidy <- tbl_graph(nodes = tax, edges = links, directed =TRUE)

routes_tidy %>% 
  activate(edges) %>% 
  arrange(desc(lineage))

ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()

#sankeyNetwork(Links = edges_d3, Nodes = nodes_d3, Source = "from", Target = "to", 
#              NodeID = "label", Value = "weight", fontSize = 16, unit = "Letter(s)")


# 5. Arbol filogenetico colorado por dos colores, barcodes parentales y absorbidos.



quit(save = 'no')


# con estos resultados cuentas los linajes completos 
.parent_id <- rownames(otu_map[otu_map$curated == "parent",])
.merged_id <- rownames(otu_map[otu_map$curated == "merged",])
parent.tax <- tax[rownames(tax) %in% .parent_id, ]
merged.tax <- tax[rownames(tax) %in% .merged_id, ]

parent.tax[complete.cases(parent.tax[c(1:5)]), ]
merged.tax[complete.cases(merged.tax[c(1:5)]), ]


rank <- 1

datavis <- rbind(data.frame(table(parent.tax[,rank]), curated = 'parent'),
      data.frame(table(merged.tax[,rank]), curated = 'merged'))


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
taxa1 <- as.matrix(tax)
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
