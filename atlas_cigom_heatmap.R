# Set colors
# x7,x5 -> misma estacion
# x4,x6 -> misma temporada

# TOMAR VALORES DE ABUNDANCIA (NO DOMINANTES) !!!
# la busqueda de los abundantes por estacion.
# Y FRECUENCIA DE DISTRIBUCION para
# 
# seleccion del top n de mas abundantes por crucero, considerar familias de interes comercial

# tamano: taxones con grupos dominantes
# 300 dpi, formato tif o png , formato tabloide doble carta,

color_Crucero <- c('X04'='#8c510a','X05'= '#bf812d', 'X06'='#dfc27d', 'X07'='red')

color_strategy <- c('OTU'='#253494', 'ASV'='#636363')
color_marker <- c('18S'='#253494', 'COI'='#636363')
Ranks1 <- c("Kingdom", "Phylum", "Class", "Order", 
           "Family", "Genus", "Species")

Ranks2 <- c("Root","Domain","Kingdom", "Phylum", "Class", "Order", 
           "Family", "Genus", "Species")
Ranks <- c("Kingdom", "Phylum", "Class", "Order", 
           "Family", "Genus", "Species")

# Set paths ----
path1 <- '/Users/cigom/metagenomics/MG_18S/multirun_xiximis/downstream_X4X5X6X7/'

p2 <- '/Users/cigom/metagenomics/MG_COI/multirun_20200507_COI/'

out_path <- '~/metagenomics/'

# functions ----

url <- 'https://raw.githubusercontent.com/RJEGR/metagenomics/master/readtx.R'

source(url)
read_rdp2 <- function(file, header = T) {
  
  tax <- read_rdp(file, header = header)
  
  tax <- tax[-8]
  
  colnames(tax) <- c(Ranks, "SL")
  
  tax[is.na(tax)] <- 'Undetermined_R'
  
  tax <-data.frame(row.names = rownames(tax),
                   mutate_all(data.frame(tax), funs(str_replace_all(., c(" "="_")))))
  return(tax)
}


# Count data ----
# 18S
c1 <- dir(path = path1, 
          pattern = '*list.shared', 
          full.names = T)

ab1 <- read.table(c1)

# COI
c2 <- dir(path = p2, 
          pattern = 'multirun_ASVs_count.table', 
          full.names = T)

ab2 <- read.table(c2)

# metadata ----
# 18S
m1 <- data.frame(
  row.names = names(ab1),
  Crucero = sapply(strsplit(names(ab1), "[_]"), `[`, 2),
  Station =  sapply(strsplit(names(ab1), "[_]"), `[`, 3),
  Strategy = '18S')

# C0I

m2 <- data.frame(
  row.names = names(ab2),
  Crucero = sapply(strsplit(names(ab2), "[.]"), `[`, 1),
  Station =  sapply(strsplit(names(ab2), "[.]"), `[`, 2),
  Strategy = 'COI')

# taxonomy ----

# 18S
t1 <- dir(path = path1, 
          pattern = '*cons.taxonomy', 
          full.names = T)
dim(tax1 <- read_rdp2(t1, header = T))

# COI
t2 <- dir(path = paste0(p2, 'Clasification'), 
          pattern = 'taxall.tsv', 
          full.names = T)

dim(tax2 <- read.delim(t2, header = F))

boots <- tax2[c(1,5,8,11,14,17,20,23,26,29)]
tax2 <- tax2[c(1,3,6,9,12,15,18,21,24,27)]

names(tax2)[2:10] <- Ranks2
names(boots)[2:10] <- Ranks2

# plot(density(rowSums(boots) / ncol(boots)))
# abline(v=median(rowSums(boots) / ncol(boots)), col = 'red')
# abline(v=mean(rowSums(boots) / ncol(boots)), col = 'blue')

# subset boots

#tax2

boots %>%
  filter(Family > 0.8) -> boots

filter_boots <- boots$V1
tax2[tax2$V1 %in% filter_boots, -1] -> tax2
ab2[rownames(ab2) %in% filter_boots, ] -> ab2

rownames(tax2) <- filter_boots
# compare reads ----

library(tidyverse)

sam1 <- names(ab1)
sam2 <- names(ab2)

# mutate_at(vars(sam),
#           function(x) {x / sum(x) * 100 }) %>%
ab1 %>%
  pivot_longer(cols = all_of(sam1), names_to = 'id_sample', 
               values_to = 'ab') %>%
  filter(ab > 1) %>%
  mutate(Marker = '18S',
         cruise = sapply(strsplit(id_sample, "[_]"), `[`, 2)) -> ab1_longer

ab2 %>%
  pivot_longer(cols = all_of(sam2), names_to = 'id_sample', 
               values_to = 'ab') %>%
  filter(ab > 1) %>%
  mutate(Marker = 'COI', 
         cruise = sapply(strsplit(id_sample, "[.]"), `[`, 1)) -> ab2_longer

lib_zise <- rbind(ab1_longer, ab2_longer)


lib_zise %>%
  mutate(ab_log2 = log2(ab+1)) %>%
  filter(ab_log2 > 10) %>%
  ggplot(aes(x=cruise, y= ab_log2, color = Marker)) +
  geom_boxplot(notch = TRUE) +
  geom_jitter(shape=19, size=1.2, alpha = 3/5,
              position=position_jitterdodge(dodge.width=0.7, 
                                            jitter.width=0.2)) +
  scale_color_manual(name=NULL,
                     values=color_marker) +
  stat_summary(
    aes(color = Marker),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4
  ) +
  theme_classic()

# 

drawr1 <- draw_rank(tax1,ab1) # rank = names(tax_table)[1]
drawr2 <- draw_rank(tax2,ab2, rank = names(tax2)[3]) # rank = names(tax_table)[1]

# Make phyloseq (ready to filter) ----

library(phyloseq)

identical(names(ab1),rownames(m1))
identical(rownames(ab1), rownames(tax1))

physeq <- phyloseq(otu_table(ab1, taxa_are_rows = TRUE),
                   tax_table(as(tax1, 'matrix')), 
                    sample_data(m1))

# COI

physeq2 <- phyloseq(otu_table(ab2, taxa_are_rows = TRUE),
                   tax_table(as(tax2, 'matrix')), 
                   sample_data(m2))

# Removing any abundance of zero and

#cleanClass <- c('Arachnida', 'Mammalia', 'Insecta','Reptilia', 'Aves')

cleanTax <- c('Undetermined_R', 'Insecta')
physeq1 %>% 
  subset_taxa(Kingdom == "Animalia") %>%
  subset_taxa(!Class %in%  cleanTax) %>%
  prune_taxa(taxa_sums(.) > 0,.) -> physeq1



# COI
cleanTax <- c('Arachnida', 'Mammalia', 'Insecta', 'Undetermined_R')
physeq2 %>% 
  subset_taxa(Kingdom == "Animals") %>%
  subset_taxa(!Class %in%  cleanTax) %>%
  prune_taxa(taxa_sums(.) > 0,.) -> physeq2


# summarize (hay algo raro aqui)
options(stringsAsFactors = FALSE)

physeq1 %>% reads_n_features(ktone = 1, Crucero == 'X04') -> x4v9
physeq1 %>% reads_n_features(ktone = 1, Crucero == 'X05') -> x5v9
physeq1 %>% reads_n_features(ktone = 1, Crucero == 'X06') -> x6v9
physeq1 %>% reads_n_features(ktone = 1, Crucero == 'X07') -> x7v9

physeq2 %>% reads_n_features(ktone = 1, Crucero == 'X04') -> x4c
physeq2 %>% reads_n_features(ktone = 1, Crucero == 'X05') -> x5c
physeq2 %>% reads_n_features(ktone = 1, Crucero == 'X06') -> x6c
physeq2 %>% reads_n_features(ktone = 1, Crucero == 'X07') -> x7c

v9Sum <- rbind(x4v9,x5v9,x6v9,x7v9)
coiSum <- rbind(x4c,x5c,x6c,x7c)
reads_n_features <- rbind(data.frame(v9Sum, Marker = '18S'), 
               data.frame(coiSum, Marker = 'COI'))

ggplot(reads_n_features, aes(x = Rank, 
                   y = ngroup, 
                   shape = Crucero, 
                   color = pct_r, 
              group = Crucero)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_path() +
  scale_color_gradient(name = 'Percent of\nSequences', low="blue", high="red") + 
  labs(caption = 'Features determined as k-tones were removed',
       x= 'Taxonomic Resolution',
       y= "Number of taxonomic groups") + 
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ Marker, scales = 'free_y')

reads_n_features %>% 
  filter(Rank == 'Family') %>%
  select(-pct_r, -pct_a)

# ABUNDANCE VS FREQUENCY


# count the number of families
taxrank <- 'Family'
top <- 50 # 50
family <- tax_glom(physeq1, taxrank = taxrank, NArm = FALSE)

otu_table <- as(otu_table(family), 'matrix') %>% as_tibble()
tax_table <- as(tax_table(family), 'matrix') %>% as_tibble()

row.names <- tax_table %>% select(all_of(taxrank))
dat <- data.frame(otu_table, 
                  row.names = row.names$Family)

total <- sum(rowSums(dat))
tax_sum <- round((rowSums(dat) / total) * 100, 4)

# boxplot ----
c <- sapply(strsplit(sam, "[_]"), `[`, 2)
s <- sapply(strsplit(sam, "[_]"), `[`, 3)
t <- substr(s, 1,1)

library(ggridges)

dat %>%
  as_tibble(rownames = 'Family') %>%
  inner_join(tax_table) %>%
  filter(Family %in% names(selection)) %>%
  # mutate_at(vars(sam), function(x) {x / sum(x) * 100 }) %>%
  pivot_longer(cols = sam) %>%
  mutate(c = sapply(strsplit(name, "_"), `[`, 2),
         s = sapply(strsplit(name, "_"), `[`, 3),
         logValue = log2(value + 0.5)) %>%
  ggplot(aes(x = Family, y = logValue), 
         fill = factor(s)) +
  geom_jitter(position=position_jitter(width=0.3, height=0.2), aes(colour=factor(s)), alpha=0.9) +
  geom_boxplot() + facet_grid(c~.) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# select top based on un-centralized samples or top of abundance fam ----
library(superheat)

dat %>% 
  as.matrix() %>% 
  superheat(col.dendrogram = T,
            row.dendrogram = T,
            clustering.method = 'hierarchical',
            dist.method = 'euclidean',
            print.plot = F ) -> sh

select_top <- head(sh$order.rows, top)

select_top <- dat[select_top,]

yr_plot_clus <- tax_sum[names(tax_sum) %in% rownames(select_top) ]
yr_plot_max <- head(tax_sum, 50)

# compare the top-abundance method:

sum(names(yr_plot_clus) %in% names(yr_plot_max))

# plot(ecdf(yr_plot_clus))
# plot(ecdf(yr_plot_max))

sam <- names(dat)

selection <- yr_plot_clus
# selection <- yr_plot_max

dat %>% 
  filter(rownames(.) %in% names(selection)) %>%
  mutate_at(vars(all_of(sam)),
          function(x) {x / sum(x) * 100 }) -> select_top_pct 


rownames(select_top_pct) <- names(selection)

superheat(ftable,
          col.dendrogram = T,
          clustering.method = 'hierarchical',
          linkage.method = 'complete',
          dist.method = 'euclidean',
          pretty.order.rows = F,
          #yr = selection,
          yr.plot.type = 'bar',
          left.label.text.size = 2.0,
          left.label.col = 'white',
          bottom.label = 'variable',
          bottom.label.text.size = 2.0,
          bottom.label.text.angle = 75,
          grid.hline = F,
          grid.vline = F,
          scale = T)

# choose top of families
topFamilies <- names(selection)

family %>% 
  subset_taxa(Family %in%  topFamilies) %>%
  prune_taxa(taxa_sums(.) > 0,.) -> phyFamilies


# This is the actual hierarchical clustering call, specifying average-link clustering
# distance methods ----


physeq.rr <- transform_sample_counts(phyFamilies, function(x) 1E6 * x/sum(x))

dist = "bray"

ord_meths = c("DCA","CCA", "RDA","NMDS", "MDS", "PCoA")

plist = plyr::llply(as.list(ord_meths), function(i, phyFamilies, dist){
  ordi = ordinate(phyFamilies, method=i, distance=dist)
  plot_ordination(phyFamilies, ordi, "samples", color="Crucero")
}, physeq.rr, dist)

names(plist) <- ord_meths

pdataframe = plyr::ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Eje_1", "Eje_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

library(ggplot2)
p = ggplot(pdataframe, aes(Eje_1, Eje_2, color=Crucero, fill=Crucero))
# p = p + geom_point(size=4)  + geom_polygon()
p = p + geom_label(aes(label = Station), fill = "white", alpha = 0.9)
p = p + facet_wrap(~method, scales="free")
#p = p + facet_grid(Crucero~method, scales="free")
p = p + scale_fill_manual(values=color_Crucero)
p = p + scale_color_manual(values=color_Crucero)
p = p + theme(legend.position = "none") +
  labs(x = "Eje 1", y = "Eje 2")

ggsave("Ordinations_topfamilies.png", p, path = p1, width = 10,height = 7)

# Networks ----
ig <- make_network(phyFamilies, max.dist=0.3, dist.fun="euclidean")
net <- plot_network(ig, phyFamilies, color="Crucero", label="Station", line_weight=0.4, alpha = 0.4) + scale_color_manual(values=color_Crucero) 
ggsave("network_topfamilies.png", net, path = p1)
# dendogram ----
library(dendextend)
library(gplots)

phyFamilies %>%
  transform_sample_counts(function(x) 1E2 * x/sum(x)) %>% 
  distance(method="euclidean", type="samples") %>%
  hclust(method="complete") %>% 
  as.dendrogram() -> dend

samples <- sapply(strsplit(labels(dend), "[_]"), `[`, 3)
cruise <- sapply(strsplit(labels(dend), "[_]"), `[`, 2)
cc <- color_Crucero[match(cruise, names(color_Crucero))]
#labels_colors(dend) <- cc
labels(dend) <- samples # paste(samples, cruise, sep = '_')

dend %>% 
  set("branches_k_color", k=3) %>%
  #rect.dendrogram(k=3, border = 8, lty = 5, lwd = 2) %>%
  set("leaves_pch", 19)  %>% 
  set("leaves_cex", 0.7) %>% 
  set("leaves_col", cc) %>% 
  #set("leaves_pch", c(17, 18, 19, 21)) %>%
  set("labels_cex", 0.7) %>%
  plot(xlab = "Samples", ylab = 'Heigh (Euclidean distance)', horiz = F)

dend %>% as.ggdend() -> ggd1

savep <- ggplot(ggd1, horiz = F, theme = theme_bw(), offset_labels = -0.01, labels_cex = 0.6) 

ggsave("dendogram_topfamilies.png", savep, path = p1)


# agglomerate rank of interest ----

agglom_rank <- function(physeq, taxrank, filter = 0.01, NArm = FALSE) {
  require(data.table)
  
  options(stringsAsFactors = F)
  
  ra <- sample_sums(physeq) > 100
  if(ra) {
    physeq <- transform_sample_counts(physeq, function(x) (x / sum (x) ) * 100)
    Rank <- tax_glom(physeq, taxrank = taxrank, NArm=NArm)
    
  } else {
    Rank <- tax_glom(physeq, taxrank = taxrank, NArm=NArm)
  }
  
  #taxaOrder = names(sort(taxa_sums(Rank), decreasing = TRUE))
  
  mdf <- psmelt(Rank)
  
  mdf <- mdf %>% filter(Abundance > 0)
  
  setDT(mdf)
  
  pos <- which(names(mdf) %in% taxrank)
  names(mdf)[pos] <- 'Group'
  
  
  mdf[, Group := ifelse(is.na(Group), 'Others', Group)]
  mdf[(Abundance < fra), Group := "Others"]
  
  return(mdf)
}

taxrank <- 'Family'
fra <- 0

# mdf <- agglom_rank(physeq1, taxrank, fra)
mdf <- phyFamilies %>%
  transform_sample_counts(function(x) 1E2 * x/sum(x)) %>%
  filter_taxa(function(x) sum(x) >= fra, TRUE) %>%
  psmelt()

dist <- distance(phyFamilies, method = "jaccard")
sampleOrder <- heatmap(as.matrix(dist))$colInd
sampleOrder <- colnames(otu_table(phyFamilies))[sampleOrder]


taxaOrder = names(sort(taxa_sums(phyFamilies), decreasing = TRUE))

mdf[mdf$Abundance == 0, 'Abundance'] <- NA

mdf$OTU <- factor(mdf$OTU, levels = taxaOrder)
mdf$Sample <- factor(mdf$Sample, levels = sampleOrder)

library(RColorBrewer)
colourCount <- length(unique(mdf$Class)) 
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
fvalues = getPalette(colourCount)

mdf %>% 
  #filter(Crucero %in% c('X04', 'X06') ) %>%
  ggplot(aes(x=Station, 
                   y=Family, 
                   color=Class)) +
  geom_point(aes(size=Abundance))+ 
  facet_grid(.~Crucero, scales="free_x") +
  scale_color_manual(values = fvalues) +
  theme_classic() 
ggsave("heatmap1_topfamilies.png", plot1, path = p1, width = 10,height = 14)

plot2 <- mdf %>% 
  filter(Crucero %in% c('X05', 'X07') ) %>%
  ggplot(aes(x=Station, 
                    y=Family, 
                    color=Class)) +
  geom_point(aes(size=Abundance))+ facet_grid(.~Crucero, scales="free_x") +
  scale_color_manual(values = fvalues) +
  theme_classic() 

ggsave("heatmap2_topfamilies.png", plot2, path = p1, width = 10,height = 14)

ggplot(mdf, aes(x = log2(Abundance+0.5), 
                y = Family, fill = Crucero)) +
  geom_density_ridges() + facet_grid(~Crucero)

# heatmap ----
#selection <- yr_plot_clus
selection <- yr_plot_max
topFamilies <- names(selection)
set <- 'X04'

f1 <- tax_glom(physeq1, taxrank = taxrank, NArm = FALSE)
f2 <- tax_glom(physeq2, taxrank = taxrank, NArm = FALSE)

f1 %>% 
  subset_taxa(Family %in%  topFamilies) %>%
  prune_taxa(taxa_sums(.) > 0,.) -> f1_top

f2 %>% 
  subset_taxa(Family %in%  topFamilies) %>%
  prune_taxa(taxa_sums(.) > 0,.) -> f2_top

pheat1 <- subset_samples(f1_top, Crucero == 'X07') %>%
  prune_taxa(taxa_sums(.) > 0, .)


pheat1 %>%
  transform_sample_counts(function(x) 1E2 * x/sum(x)) -> heatFam

dist <- distance(heatFam, method = "euclidean")
sampleOrder <- heatmap(as.matrix(dist))$colInd
sampleOrder <- colnames(otu_table(heatFam))[sampleOrder]

taxaOrder = names(sort(taxa_sums(heatFam), decreasing = TRUE))

#
dist %>%
  hclust(method="complete") %>% 
  as.dendrogram() -> dend
#
plot(dend)

otu_table(heatFam)[otu_table(heatFam) == 0] <- NA

#taxaOrder = rownames(sort(unique(tax_table(Family)[,'Orden']), decreasing=TRUE))
colourCount <- 8 
getPalette = colorRampPalette(brewer.pal(colourCount, "YlGnBu"))
fvalues = getPalette(colourCount)

h <- NULL
h <- plot_heatmap(heatFam, 
                  sample.label = "Station", 
                  sample.order = sampleOrder,
                  taxa.label = taxrank,
                  taxa.order = taxaOrder,
                  na.value = "white", trans = NULL,
                  low = fvalues[1:4],
                  high = fvalues[4:8]
) +
  labs(fill = "Abundancia\nRelativa (%)") +
  theme(legend.position = "left") +
  geom_tile(color = "grey")

h <- h + guides(fill = guide_colorbar(barheight = unit(7, "cm"),  ticks.colour = "black", frame.colour = "black"))

plot <- NULL

plot <- h + facet_grid(Phylum+Class ~ ., 
                       scales = "free", space = "free" 
                       #                             switch = "y"
) + 
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1),
    strip.text.y = element_text(
      angle = 0, 
      size = 5),
    strip.background = element_rect(colour = "black", 
                                    fill = "transparent",
                                    size = 0.1),
    panel.spacing = unit(0.001, "lines")
  ) +
  labs(fill = "Abundancia\nRelativa (%)")


ggsave("heatmap1_topfamilies.png", plot, path = p1, width = 10,height = 14)


# simper ----

library("FactoMineR")
library('FactoInvestigate')
library("factoextra")
library("corrplot")

ss <- f1_top
otu_table <- as(otu_table(ss), 'matrix') %>% t()
tax_table <- as(tax_table(ss), 'matrix') %>% as_tibble() %>% select(all_of(taxrank))

colnames(otu_table) <- tax_table$Family

# metadat
sam <- rownames(otu_table)
# c <- sapply(strsplit(sam, "[_]"), `[`, 2)
# s <- sapply(strsplit(sam, "[_]"), `[`, 3)
t <- substr(s, 1,1)

c <- sapply(strsplit(sam, "[.]"), `[`, 1)
s <- sapply(strsplit(sam, "[.]"), `[`, 2)

dat <- data.frame(SampleID = sam, otu_table, ship = s, cruise = c, transect = t)

#datt <- t(dat)
rownames(dat) <- paste(s,c, sep = "_")

pca <- dat %>%
  select(-ship, -cruise, -transect, -SampleID) %>%
  PCA(., graph = FALSE)

# cluster

set.seed(080620)
var <- get_pca_var(pca)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)

fviz_pca(pca, geom = c("point", "text"), col.ind = c,  palette = "jco" ,label = "ind", geom.var = "point", repel = T)

fviz_pca(pca, geom = c("point", "text"), col.ind = c,  palette = "jco" ,label = "var", geom.var = "text", repel = T) # col.var = grp

fviz_eig(pca, addlabels = TRUE, ylim = c(0, 20))

var <- get_pca_var(pca)
corrplot(var$cos2, is.corr=FALSE)
corrplot(var$contrib, is.corr=FALSE)
# fviz_pca_var(pca, alpha.var = "cos2")

fviz_pca_var(pca, 
             col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE,
             geom = "text",
             legend.title = list(col = "quality of\nrepresentation"))

#

set.seed(080620)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)

# Color variables by groups
fviz_pca_var(pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")


ind.p <- fviz_pca_var(pca, 
                      geom = "text", 
                      col.ind = c )


ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis",
              subtitle = "Main taxa reported in 18S",
              caption = "",
              ticks = TRUE,
              xlab = "PC1", ylab = "PC2",
              legend.title = "Ship", 
              legend.position = "top",
              ggtheme = theme_gray(), palette = "jco" )


fviz_pca_var(pca, 
                col.ind = c, 
                palette = "jco", 
                addEllipses = F, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Ship",
                geom.var = c("text"))


fviz_pca_biplot(pca, 
                # Individuals
                geom.ind = "text",
                geom.var = 'text',
                fill.ind = dat$cruise, 
                #col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                #repel = T,
                addEllipses = F,
                # Variables
                #alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdYlBu",
                legend.title = list(fill = "Ship", color = "Contrib",
                                    alpha = "Contrib"))

#

res.hcpc = HCPC(pca, nb.clust = -1, graph = FALSE)

fviz_dend(res.hcpc, 
          cex = 0.7,  
          palette = "jco", 
          rect = TRUE, rect_fill = TRUE, 
          rect_border = "jco",  
          labels_track_height = 0.8 )


fviz_cluster(res.hcpc,
             repel = TRUE, geom = 'point',
             show.clust.cent = TRUE, 
             palette = "jco",
             ggtheme = theme_minimal(),
             main = "Factor map"
)

# go to the lab ----

