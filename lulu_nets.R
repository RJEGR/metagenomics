#!/usr/bin/env Rscript
rm(list=ls())
setwd(dir="/Users/cigom/metagenomics/LULU")
# for i in $(seq 0.010 0.002 0.052); do Rscript --vanilla lulu_nets.R $i; done
# :::::::::::::::::: loading args

args = commandArgs()
cutoff = args[7]

#cutoff = 0.03

# ::: load libs

library(igraph)
library(ggraph)
library(tidygraph)
library(tidyverse)

# :::::::::::::::::: loading data after processing lulu_cutoff.R
# 
temp = list.files(path = "./otu_map", pattern="*.otu_map.csv", full.names = TRUE)
myfiles = lapply(temp, read.csv, sep = " ", header=TRUE, stringsAsFactors=FALSE)
myfiles <- lapply(myfiles, rownames_to_column, var="id")
data <- data.frame(do.call(rbind, myfiles))

# ::::::::::::::::: loading constaxonomy and shared through ls
temp <- NULL
temp = list.files(path = "./shared.tax", pattern="shared.taxonomy.csv", full.names = TRUE)
myfiles2 = lapply(temp, read.csv, sep = " ", header=TRUE, stringsAsFactors=FALSE)
myfiles2 <- lapply(myfiles2, rownames_to_column, var="id")
data2 <- data.frame(do.call(rbind, myfiles2))
data2[is.na(data2$Especie) , 'Especie'] <- "Undetermined"


# :::::::::::::::::::: loading curated dataset of otus 
temp <- NULL
temp = list.files(path = "./curated_table", pattern="curated_table.csv", full.names = TRUE)
myfiles3 = lapply(temp, read.csv, sep = " ", header=TRUE, stringsAsFactors=FALSE)
myfiles3 <- lapply(myfiles3, rownames_to_column, var="id")
data3 <- data.frame(do.call(rbind, myfiles3))

# ===================
# ===================
# ===================

# buscamos el otu en el que se representan el mayor numero de merged

# ::: select the highest merged otu
data %>%
  as.tibble() %>%
  filter(!!cutoff == cutoff) %>%
  select(parent_id) %>%
  group_by(parent_id) %>%
  summarize(n=n()) %>%
  arrange(desc(n)) %>% # sample_n(1)
  slice(1L) %>% # select the position 1 as the highest otu
  inner_join(data, by = "parent_id") %>%
  filter(!!cutoff == cutoff) -> highest_otu



# NETWOKR VISUALIZATION

highest_otu %>%
    select(parent_id) %>%
    unique() -> otu_name

highest_otu %>%
    select(id) %>%
    unique() -> otu_list

vsearch <- read.table(paste0("/Users/cigom/metagenomics/LULU/vsearch/vsearch.", cutoff, ".list.txt"), header=FALSE)
links <- vsearch
colnames(links) <- c("id","node2","dist")

links %>%
   as.tibble() %>% 
   inner_join(otu_list, by = "id") -> merged_net

# parent_net <- links[links$id == paste(otu_name),]

net <- data.frame(from = merged_net$id, 
                    to = merged_net$node2, 
                    dist = merged_net$dist)

# merged_net <- filter(data, curated == "merged" ,cutoff == c("0.03", "0.05"))

#net <- data.frame(from = merged_net$id, 
#                    to = merged_net$parent_id, 
#                    value = merged_net$spread, 
#                    cutoff = merged_net$cutoff)

# Number of connection per OTU
c( as.character(net$from), as.character(net$to)) %>%
  as.tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> nodes
colnames(nodes) <- c("id", "n")

# :::::::::::::::::::::
# Add, Classification 
# ::::::::::::::::::::

data2 %>%
    as.tibble() %>%
    filter(!!cutoff == cutoff) %>%
    inner_join(nodes, by = "id") %>%
    select(id, n, Orden, Familia, Genero, Especie, cutoff) %>%
    arrange(id) -> nodes



nodes$merged <- FALSE
nodes[nodes$id %in% otu_list$id, 'merged' ] <- TRUE
# Rename NAs
nodes[is.na(nodes$Especie) ,'Especie'] <- "Undetermined"
table(nodes$Especie)


write.table(nodes, file = paste0('nodes.',cutoff,'.csv'),
                sep = " ", quote = FALSE)


quit(save = "no")

data2[is.na(data2$Familia) , 'Familia'] <- "Undetermined"

data2 %>%
  as.tibble() %>%
  filter(Familia == "Undetermined") %>%
  select(cutoff, Familia) %>%
  group_by(cutoff) %>%
  summarize(Family=n()) -> Undetermined_f

data2 %>%
  as.tibble() %>%
  filter(Especie == "Undetermined") %>%
  select(cutoff, Especie) %>%
  group_by(cutoff) %>%
  summarize(Species=n()) -> Undetermined_e

plot_und <- cbind(Undetermined_f, Undetermined_e[,2])



9,790 / 16,936

# ::::::::::
# re-open data
# and start wrapping results
# ===========
temp <- NULL
myfiles <- NULL
temp = list.files(path = "./nodes_by_cutoff/", pattern="nodes.*.csv", full.names = TRUE)
myfiles = lapply(temp, read.csv, sep = " ", header=TRUE, stringsAsFactors=FALSE)
data4 <- data.frame(do.call(rbind, myfiles))

library(ggplot2)
#
ggplot(data4[data4$merged == "TRUE", ], aes(x=id, y=n, group=Especie)) +  
    #geom_line(colour="grey") +
    #geom_line(aes(linetype=Especie)) +
    geom_point(aes(color=Especie, shape=Especie), size=1.5, alpha=0.3) +
    # geom_point(aes(size=n, color=merged, shape=Especie)) + 
    scale_color_brewer(palette = "Set1")  +
    theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background =element_blank() ) +
    facet_wrap(~ cutoff, scales="free_y")

# ======

ggplot(data4[data4$merged == "FALSE", ], aes(log10(n), ..count.., colour=Especie, fill=Especie)) +
    #geom_histogram(alpha=0.35, position = "stack", bins = 500)  +
    geom_density(alpha=0.35, position = "stack") +
    scale_fill_brewer(palette = "Set1" ) +
    scale_color_brewer(palette ="Set1" ) +
    facet_wrap(~ cutoff, scales="free_y") + theme_minimal()

# ========= NEXT

data4 %>%
    as.tibble() %>%
    filter(cutoff == 0.03) %>%
    select(id, merged) -> data4_subset


data2 %>%
    as.tibble() %>%
    filter(cutoff == 0.03) %>%
    inner_join(data4_subset, by = "id") %>%
    select(id, Especie, paste(colnames(data2))[c(2:48)], merged) %>%
    arrange(Especie) -> parent_species

parent_species_melt <- reshape2::melt(parent_species)

library(ggpubr)

ggbarplot(parent_species_melt[parent_species_melt$merged == "TRUE", ], x = "variable", y = "value",
          fill = "Especie", group = "Especie",
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          x.text.angle = 90,           # Rotate vertically x axis texts
          facet.by = "Especie")


ggplot(parent_species_melt[parent_species_melt$merged == "TRUE", ], aes(log10(value), ..count.., colour=Especie, fill=Especie)) +
    #geom_histogram(alpha=0.35, position = "stack", bins = 500)  +
    geom_density(alpha=0.35, position = "stack") +
    scale_fill_brewer(palette = "Set1" ) +
    scale_color_brewer(palette ="Set1" ) +
    facet_wrap(~ variable, scales="free_y") + theme_minimal()


# el heatmap no es una buena idea debido a que hay muchos valores de cero
# =========

matrix <- data3[, -c(1, ncol(data3))]
#library(vegan)
specnumber <- aggregate(x = matrix, by = list(data3$cutoff), FUN = c("specnumber")) # observed number of species (Beta-diversity)
diversity <- aggregate(x = matrix, by = list(data3$cutoff), FUN = c("diversity")) #shannon

rich <- t(specnumber[,-c(1)])
colnames(rich) <- specnumber[,c(1)]

chisq <- chisq.test(rich)
library(corrplot)
corrplot(chisq$residuals, is.cor = FALSE)


# ===
# merged <- data.frame(shared_taxa[!shared_taxa$id %in% parent$id, ])
# ===

# Create a graph object with igraph
mygraph <- graph_from_data_frame(net, vertices = nodes, direct = TRUE)

netplot <- igraph::simplify(mygraph, remove.multiple=  TRUE, remove.loops = TRUE)
graph <- as_tbl_graph(netplot) # %>% 
    #mutate(Connections = centrality_degree(mode = 'in'))


png(filename = paste0('network.',cutoff,'.png'),  width = 480, height = 380, units = "px",  res = 100)

ggraph(graph, layout="igraph", algorithm="kk") + 
    geom_edge_fan(aes(alpha = ..index..), show.legend = FALSE) + 
    geom_node_point(aes(shape = merged, color = Especie , size=n)) + 
    #facet_edges(~ merged) + 
    scale_color_brewer(palette = "Set1") +
    theme_graph(foreground = 'steelblue', fg_text_colour = 'white') #+
    #facet_nodes(~ merged)

dev.off()

quit(save = "no")

# ===============
# : 
# ===============

# Make the graph

#ggraph(mygraph, layout="igraph", algorithm = 'kk') + 
#  geom_edge_link(edge_colour="black", edge_alpha=0.2, edge_width=0.3) +
#  geom_node_point(aes(size=n, alpha=n, color = merged)) +
#  scale_fill_viridis(palette = "Paired") +
#  theme_void() +
#  theme(
#    legend.position="none",
#    plot.margin=unit(rep(1,4), "cm")) +
#    facet_nodes(~ Especie)

#  ====================


# further analysis
# extraemos estos otus del objeto data2
data2 %>%
  as.tibble() %>%
  filter(cutoff == paste(cutoff)) %>%
  inner_join(highest_otu, by = "id") -> data2_subst




# preparamos matrix para graficar 
from <- ncol(data2_subst) - 14
to <- ncol(data2_subst)
matrix <- data2_subst[, -c(1, from:to)]
# correlacionamos datos


#heatmap(cor(matrix))
#M <- cor(matrix, method = "pearson")
p.mat <- cor.mtest(matrix)
#library(corrplot)
#corrplot(M, type="upper", order="hclust", col=c("black", "white"),
#         bg="lightblue",
#         tl.col="black", tl.srt=45,
#         p.mat = p.mat, sig.level = 0.05, insig = "blank")

# ===
# maxs <- apply(matrix, 2, max) # max values over the columns
# mins <- apply(matrix, 2, min)
# library(scales)
# normalized <- as.data.frame(scale(matrix, center = mins, scale = maxs-mins))
# datanorm <- cbind(matrix, data2_subst[, c(1, (from+1):(from+7))])


# === melt data
library(reshape2)

dataplot <- melt(datanorm, 
                variable.name = "Estacion",
                value.name = "Abundancia")
# === plot
library(ggplot2)
ggplot(dataplot, aes(Estacion, log10(Abundancia), group=Estacion) ) +
  geom_line() +
  geom_point(aes(color=Especie)) +
  coord_flip()
  
# not good!
library(ggridges)
ggplot(dataplot, aes(x = Abundancia, y = Estacion, fill = Especie)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")


# == boxplot
ggplot(dataplot, aes(x=Estacion, y=Abundancia)) + 
geom_boxplot() + 
stat_summary(fun.y = mean, geom="point",colour="black", size=1.5) +
theme(
        axis.text.x = element_text(angle=90, hjust=1))
# == density plot

ggplot(dataplot, aes(log10(Abundancia), colour=Especie, fill=Especie)) + 
    geom_density() + facet_wrap(~Estacion)

# Chi-square test basics: 
# Chi-square test examines whether rows and columns of a contingency table are statistically significantly associated.

chisq <- chisq.test(t(matrix))
library(corrplot)
corrplot(chisq$residuals, is.cor = FALSE)
# Contibution in percentage (%)
# The contribution (in %) of a given cell to the 
# total Chi-square score is calculated as follow:
contrib <- 100*chisq$residuals^2/chisq$statistic
autoplot(as.matrix(chisq$residuals))

pca <- prcomp(contrib, scale=TRUE)
pca.var <- pca$sdev^2

library("factoextra")
fviz_pca_biplot(pca, col.ind = data2_subst$Especie,
                palette = "jco", geom = "point")

# ===


# ===================
# ===================
# ===================
# ===================













# :::::::::::::::::::::::
# https://www.udemy.com/data-science-and-machine-learning-bootcamp-with-r/?gclid=CjwKCAjwmJbeBRBCEiwAAY4VVZjhYmHbrpq8PEa2DbBPzQcn0dbAi2T-K35O2v-_BuXs9ntlS6F7RxoCWQgQAvD_BwE&k_clickid=a51505da-84f6-4331-b84c-1ad3bb47479b_408_GOOGLE_NEW-AW-PROS-PROF-Bus-DataScience-EN-ENG_._ci_821726_._sl_ENG_._vi_PROF_._sd_All_._la_EN_.__topic_%2Bdata+%2Bscience_b_281634940786_c&utm_medium=udemyads&utm_campaign=NEW-AW-PROS-PROF-Bus-DataScience-EN-ENG_._ci_821726_._sl_ENG_._vi_PROF_._sd_All_._la_EN_._&utm_term=_._pl__._pd__._ti_kwd-54109157547_._kw_%2Bdata+%2Bscience_._&utm_content=_._ag_topic_._ad_281634940786_._de_c_._dm__._lo_9073793_._&utm_source=adwords&matchtype=b
# Neural nets (matchine learling approach)
# scaling data (Normalize)
# :::::::::::::::::::::::

data3_sub <- filter(data3, cutoff == c("0.03"))

tmp <- data3_sub[, -c(1, ncol(data3_sub))]
total <- rowSums(tmp)
tmp <- cbind(total, tmp)

maxs <- apply(tmp, 2, max ) # max values over the columns
mins <- apply(tmp, 2, min)

library(scales)
# center column to the minimun values and devide values to the max values
scale.data3 <- as.data.frame(scale(tmp, center = mins, scale = maxs-mins))

#install.package('caTools')
library(caTools)
split <- sample.split(scale.data3, splitRatio = 0.7)
train <- subset (scale.data3, split == TRUE)
test <- subset (scale.data3, split == TRUE)

library(neuralnet)
n <- names(train)

# ====== b-diversity


beta <- data.frame(Samples = rownames(t(scale.data3)) , t(scale.data3))
beta2 <- plyr::ddply(beta, ~Samples, function(x) { data.frame(RICHNESS=sum(x[-1]>0))})            


f <- as.formula(paste("total ~", paste(n[!n  %in% "total"], collapse = " + ")))

nn <- neuralnet(f, data=train, hidden = c(5, 3), linear = TRUE)
plot(nn)

# :: prediction

predicted.nn.values <- compute(nn, test[,-c(1)])

true.predictions <- predicted.nn.values$net.result * (max(tmp$total)-min(tmp$total))+min(tmp$total)

# convert test data
test.r <- (test$total) * (max(tmp$total) - min(tmp$total)) + min(tmp$total)
MSE.nn <- sum( test.r - true.predictions^2) / nrow(test)
error.df <- data.frame(test.r, true.predictions)

library(ggplot2)
ggplot(error.df, aes(x=test.r, y= true.predictions)) +
    geom_point() + stat_smooth()

# repetir esta prueba con los datos de 
# alfa y beta diversidad y 
# algun otro valor a lo largo de las estaciones.