library(tidyverse)
library(igraph)


dir <- "/Users/cigom/Desktop/STRING_Blanca/"
setwd(dir)
k <- read.table("string_interactions.tsv", header=FALSE)
colnames(k) <- c("node1", "node2", 
        "node1_string_internal_id",
        "node2_string_internal_id",
        "node1_external_id",       
        "node2_external_id",
        "neighborhood_on_chromosome",
        "gene_fusion",
        "phylogenetic_cooccurrence",       
        "homology",        
        "coexpression",      
        "experimentally_determined_interaction", 
        "database_annotated",     
        "automated_textmining",    
        "combined_score")

links <- k[,c("node1", "node2","combined_score")]


#::: prepare list to compare versus metadata
list <- c(as.vector(links$node2), as.vector(links$node1))
list <-list[order(list)]
list <- list[!duplicated(list)]

# ::::::::::::::
# prepare adjacense matrix
# :::::::::::::::



# and correlate by expression size:
# cut -d"," -f1-4,6 nodes_relabel.csv > nodes_relabel.csv.in
nodes <- read.table("nodes_relabel.csv.in", sep=",", stringsAsFactors=FALSE, header=TRUE)
colnames(nodes)[c(1,6:7)] <- c("Contig","id", "Node")
nodes$Specie <- sapply(strsplit(nodes$id, "_"), "[", 2)
# nodes[389,4] <- 0
nodes <- na.omit(nodes)
nodes <-nodes[order(nodes$Node),]

# :::
# and normalize values by sample
# :::

n <- nodes[,c(2:5)]
#nodes[,c(2:5)] <- n / rowSums(n)
nodes[,c(2:5)] <- log2(n+1)

png(filename = "log2_plus1_violin.png",  
    width = 12, height = 6, units = "in",
    res = 300)
ggpubr::ggviolin(reshape::melt(log2(n+1)), 
                x = "variable", 
                y = "value", 
                fill = "variable"
                )
dev.off()
png(filename = "violin.png",  
    width = 12, height = 6, units = "in",
    res = 300)
ggpubr::ggviolin(reshape::melt(n), 
                x = "variable", 
                y = "value", 
                fill = "variable",
                #yscale = "log2"
                )
dev.off()

# label again some proteins
nodes$Node[c(21,37,45, 66, 181)] <- c("ANGEL1", "BHMT", "CALM2", "ACADSB", "KRT1")

#:::
T <- list[list %in% nodes$Node]
F <- list[!list %in% nodes$Node]
# ::: get cellular function
cf <- read.table("Genes_DEGs_cell_funct.tsv", sep="\t", stringsAsFactors=FALSE, header=TRUE)
cf <-cf[order(cf$Entry.name),]
cf <-cf [!duplicated(cf$Entry.name),]
cf$Entry.name[c(19,38,47, 67, 186)] <- c("ANGEL1", "BHMT", "CALM2", "ACADSB", "KRT1")

# ::
# merge edge labels 
# ::
x <- data.frame(nodes[nodes$Node %in% list,c("Node", "Contig","BC", "BV", "CC", "CV")], row.names=seq(1, length(list)))
x <-x[order(x$Contig),]
cf <-cf[order(cf$contig),]
y <- cf[cf$contig %in% x$Contig,]
y <- y[,3]
v <- data.frame(cbind(x, y), row.names=seq(1, length(list)))
v$Contig <- NULL

# homogenize the column y
v$y <- paste(toupper(substring(v$y, 1,1)), substring(v$y, 2), sep="")

# :: Make igraph object
net <- graph_from_data_frame(d=links, vertices=v, directed=TRUE)

# ::::::::::::::::::::::
# properties of topology
# ::::::::::::::::::::::

adj_mat<-as.matrix(get.adjacency(net))

#
pdf(file="Network_topology.pdf")
plot(degree_distribution(net),
    #log="xy",
    col="sienna4",
    xlab = "Number of connections",
    ylab = "Density")

boxplot(degree(net), col="snow2")

hist(degree(net, mode="in"),col = "tomato", 
    xlab = "Degree")

hist(degree(net, mode="out"),col = "tomato", 
    xlab = "Degree")

hist(degree(net, mode="all"),col = "tomato", 
    xlab = "Degree")

heatmap(adj_mat, Rowv=NA, Colv="Rowv")

dev.off()


net <- igraph::simplify(net, remove.multiple=  FALSE, remove.loops = TRUE)


# V(net)$label <- rownames(v)

# ::
# Generate colors based on celullar funct:
# ::

c <- levels(as.factor(V(net)$y))
# display.brewer.all()

colrs <- RColorBrewer::brewer.pal(n =length(c), name = 'BrBG')
colrs[which(c=="Non reported")] <- "gray80"

# :::
#improve plot of interest
# :::

head(igraph::as_data_frame(net, what="vertices"))

# Compute node degrees (#links) and use that to set node size:
# deg <- degree(net, mode="all")
# V(net)$size <- deg*3

# ::
# change arrow size and edge color:
# E(net)$arrow.size <- .2
# ::
# Generate edge color variable to plot the path:
new.paths <- shortest_paths(net, 1:ecount(net), output = "both")
new.paths <- ifelse( strength(net)>=4, V(net)$name, NA ) 
ecol <- rep("gray80", ecount(net))
ecol[unlist(new.paths)] <- "orange"

# E(net)$edge.color <- "gray80"

# :::
# attributes to edges
# :::

str(edge_attr(net))

# ::
# :: change the width of the edges connectivity
E(net)$width <- log(E(net)$combined_score) +1

# ::
# :: Filter the nodes-labels by number of label-interactions
V(net)$label <- ifelse( strength(net)<5, V(net)$name, NA ) 



# ::
# :: Select the the layout 
lay <- layout_nicely(net)
lay <- lay *0.8
#lay <- layout.fruchterman.reingold(net,
#            niter=1000,
#            minx=strength(net)*0.5,
#            miny=strength(net)*0.5)

png(filename = "network_1.png",  
    width = 8, height = 8, units = "in",
    res = 460)
plot(net, 
  edge.arrow.mode=0,
  #edge.color=ecol, 
  layout=lay, 
  main="String network", 
  vertex.color=colrs,
  vertex.label=V(net)$name,
  #edge.curved=.25,
  vertex.size=10,
  vertex.label.cex = .37,
  vertex.frame.color="#ffffff",
  vertex.label.color="black",
  #rescale=FALSE
    )
legend(x=-1.5, y=-1.1, paste0(c), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=2)

dev.off()

# ::
# second
# ::

pdf(file="samples_network_2_5edges_less.pdf")
mar <- c(1,1,1,1)
par(mfrow=c(2,2), mar=mar-0.2)
plot(net, edge.arrow.mode=0, 
  layout=lay, 
  main="Carrizal-Verde", 
  vertex.color=colrs,
  #vertex.label=NA,
  #edge.curved=.25,
  #vertex.size=20,
  vertex.size=V(net)$CV,
  vertex.label.cex = .27,
  vertex.frame.color="#ffffff",
  vertex.label.color="black")

# 2
plot(net, edge.arrow.mode=0, 
  layout=lay, 
  main="Carrizal-Cafe", 
  vertex.color=colrs,
  #vertex.label=NA,
  #edge.curved=.25,
  #vertex.size=20,
  vertex.size=V(net)$CC,
  vertex.label.cex = .27,
  vertex.frame.color="#ffffff",
  vertex.label.color="black")

plot(net, edge.arrow.mode=0, 
  layout=lay, 
  main="Boquita-Cafe", 
  vertex.color=colrs,
  #vertex.label=NA,
  #edge.curved=.25,
  #vertex.size=20,
  vertex.size=V(net)$BC,
  vertex.label.cex = .27,
  vertex.frame.color="#ffffff",
  vertex.label.color="black")

plot(net, edge.arrow.mode=0, 
  layout=lay, 
  main="Boquita-Verde", 
  vertex.color=colrs,
  #vertex.label=NA,
  #edge.curved=.25,
  #vertex.size=20,
  vertex.size=V(net)$BV,
  vertex.label.cex = .27,
  vertex.frame.color="#ffffff",
  vertex.label.color="black")

dev.off()



# can also remove networks by
# file:///Users/cigom/Documents/GitHub/IIIBSS/Networks/Exercises/Exercises_4_Network_Robustness.html


png(filename = "_network_1_gray.png",  
    width = 10, height = 10, units = "in",
    res = 460)
plot(net, 
  edge.arrow.mode=0,
  #edge.color=ecol, 
  layout=lay, 
  main="String network", 
  vertex.color="gray99",
  vertex.label=V(net)$name,
  #edge.curved=.25,
  vertex.size=10,
  vertex.label.cex = .37,
  vertex.frame.color="#E0E0E0",
  vertex.label.color="black",
  #rescale=FALSE
    )
dev.off()


samples <- c("CV", "CC", "BC","BV")
samples <- as.data.frame(samples)
for (i in 1:nrow(samples))  {
    print(V(net)$samples[i,1])
}

pdf(file="samples_network_2_gray.pdf")
par(mfrow=c(2,2), mar=c(1,1,1,1))
plot(net, edge.arrow.mode=0, 
  layout=lay, 
  main="Carrizal-Verde", 
  vertex.color="#E0E0E0",
  vertex.label=V(net)$name,
  vertex.size=V(net)$CV,
  vertex.label.cex = .37,
  vertex.frame.color="#B2182B",
  vertex.label.color="black")

  plot(net, edge.arrow.mode=0, 
  layout=lay, 
  main="Carrizal-Cafe", 
  vertex.color="#E0E0E0",
  vertex.label=V(net)$name,
  vertex.size=V(net)$CC,
  vertex.label.cex = .37,
  vertex.frame.color="#B2182B",
  vertex.label.color="black")

  plot(net, edge.arrow.mode=0, 
  layout=lay, 
  main="Boquita-Cafe", 
  vertex.color="#E0E0E0",
  vertex.label=V(net)$name,
  vertex.size=V(net)$BC,
  vertex.label.cex = .37,
  vertex.frame.color="#B2182B",
  vertex.label.color="black")

  plot(net, edge.arrow.mode=0, 
  layout=lay, 
  main="Boquita-Verde", 
  vertex.color="#E0E0E0",
  vertex.label=V(net)$name,
  vertex.size=V(net)$BV,
  vertex.label.cex = .37,
  vertex.frame.color="#B2182B",
  vertex.label.color="black")

  dev.off()


  # try clustering
cluster <- cluster_edge_betweenness(net)

png(filename = "cluster_egde_betweenness.png",  
    width = 10, height = 10, units = "in",
    res = 300)
plot(cluster, net,
  edge.arrow.mode=0,
  layout=lay, 
  main="Edges betweenness", 
  #vertex.color="gray99",
  vertex.label=V(net)$name,
  vertex.size=10,
  vertex.label.cex = .37,
  vertex.frame.color="#E0E0E0",
  vertex.label.color="black"
     )
dev.off()
#
# ::::::::::::::::::
#

library(tidygraph)
library(ggraph)
#  https://www.jessesadler.com/post/network-analysis-with-r/
routes_tidy <- tbl_graph(nodes = v, edges = links, directed =TRUE)
# or
# routes_igraph_tidy <- as_tbl_graph(net)
class(routes_tidy)
# ex
routes_tidy %>% 
  activate(edges) %>% 
  arrange(desc(combined_score))

# then plot
# basic topology
ggraph(routes_tidy) + geom_edge_link() + geom_node_point() + theme_graph()

#
ggraph(routes_tidy, layout = "nicely") + 
  geom_node_point() +
  geom_edge_link(aes(width = combined_score), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = Node), repel = TRUE) +
  # labs(edge_width = "Letters") +
  theme_graph()

ggraph(routes_tidy, layout = "linear") + 
  geom_edge_arc(aes(width = combined_score), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = Node)) +
  labs(edge_width = "Edges") +
  theme_graph()

