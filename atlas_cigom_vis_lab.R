
library("FactoMineR")
library('FactoInvestigate')
library("factoextra")
library("corrplot")
library(tidyverse)

color_Crucero <- c('X04'='#00AFBB','X05'= '#E7B800', 'X06'='#FC4E07', 'X07'='red')

cleanTax <- c('Arachnida', 'Mammalia', 'Insecta', 'Undetermined_R')

Level <- 'Family'

# agglomerate families
tbl_agglom <- function(tax_tbl, count_tbl, Level = names(tax_tbl)[4], 
                       fkingdom = 'Animals') {
  
  scheck <- identical(rownames(tax_tbl), rownames(count_tbl))
  
  if(!scheck) {
    sortt <- rownames(tax_tbl)
    count_tbl <- count_tbl[match(rownames(count_tbl), sortt),]
  }
  
  tbl <- cbind(tax_tbl, count_tbl)
  sam <- colnames(count_tbl)
  
  cleanClass <- c('Insecta', 'Mammalia', 'Arachnida')
  
  tbl_agglom <- tbl %>%
    filter(Kingdom == fkingdom) %>%
    filter(!Class %in%  cleanClass) %>%
    #select(Level, sam) %>%
    #drop_na(Level) %>% # to clean the redundant names (NA assignation)
    group_by_at(vars(one_of(Level))) %>%
    summarise_at(vars(sam), sum) %>%
    ungroup() %>%
    filter(Family != 'Undetermined_R')
    
  
  return(tbl_agglom)
}

 
v9_fam <- tbl_agglom(tax1, ab1, Level = Level, fkingdom = 'Animalia')
# coi_fam <- tbl_agglom(tax2, ab2, Level = Level, fkingdom = 'Animals')

total <- sum(rowSums(v9_fam[-1]))
tax_sum <- round((rowSums(v9_fam[-1]) / total) * 100, 4)

# coi_fam %>%
#   select(-Level) %>%
#   mutate(total = sum(rowSums(.))) %>% 
#   mutate(pct = round((rowSums(.) / total) * 100, 4)) %>%
#   select(pct)

# transform-data-count, how?
library(ggforce)

sam <- colnames(v9_fam)[-1]
v9_fam %>% 
  # mutate_at(vars(!all_of(Level)), function(x) {x / sum(x) * 100 }) %>%
  pivot_longer(cols = !all_of(Level), names_to = 'id_sample', 
               values_to = 'ab') %>%
  filter(ab > 1) %>% # Not!!
  mutate(cruise = sapply(strsplit(id_sample, "[_]"), `[`, 2),
         sample = sapply(strsplit(id_sample, "[_]"), `[`, 3),
         Family = reorder(Family, ab)) -> dv 
dv %>%
  ggplot(aes(x = as.numeric(Family), y = ab, fill = cruise)) +
  geom_bar(stat="identity",
           position = position_dodge(), width = 0.7,
           alpha = 0.8, color = 'black', size = 0.2) + 
  #facet_grid(cruise~.) +
  scale_fill_manual(values=color_Crucero) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  facet_zoom(x = ab > 10000) +
  scale_x_continuous(
  breaks = 1:length(levels(dv$Family)),
  label = levels(dv$Family))



# frequency
dv %>%
  filter(ab > 0) %>%
  group_by(Family, sample) %>%
  summarise(mean = mean(ab), Total = sum(ab), 
            Freq = n(), sd = sd(ab), var = var(ab)) %>%
  ungroup() %>%
  arrange(desc(Total)) -> dv_summary


dv_summary %>% 
  mutate(mean_group = ifelse(mean < 1E2, '< 1E2', 
                        ifelse(mean >= 100 & mean < 1E3, 
                               ' >= 100|< 1E3',
                               ifelse(mean >= 1E3 & mean < 2E4, '>= 1E3|< 2E4', 
                                      ifelse(mean >= 2E4 & mean < 4E4, '>= 2E4|< 4E4', 
                                             '> 4E4'))))) -> dv_summary
dv_summary %>%
  filter(Family != 'Undetermined_R') %>%
  ggplot(aes(x = Family, y = sample)) +
  theme_classic() +
  geom_point(aes(color = mean_group), #shape = factor(Freq)
             alpha = 0.7, size = 1, shape = 15) + # shape = 15
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
        axis.text.y = element_text(size = 3.7),
        legend.position = "top") +
  scale_color_brewer(type='qual', palette = 'Paired', direction = -1) +
  facet_grid(factor(Freq)~., scales = "free_y")
  # facet_zoom(x = mean > 1E3) +
  # scale_x_continuous(
  #   breaks = 1:length(levels(dv$Family)),
  #   label = levels(dv$Family))

# PCA ----
tax <- fam_top$Family
f <- fam_top %>% select(-Level)
ftable <- as(f, 'matrix') %>% t()
colnames(ftable) <- tax

# metadat
sam <- colnames(f)
# For 18S
c <- sapply(strsplit(sam, "[_]"), `[`, 2)
s <- sapply(strsplit(sam, "[_]"), `[`, 3)
t <- substr(s, 1,1)


rownames(ftable) <- paste(s,c, sep = "_")

#datt <- t(dat)


pca <- ftable %>%
  PCA(., graph = FALSE)


set.seed(080620)
var <- get_pca_var(pca)
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)

fviz_pca(pca, geom = c("point", "text"), col.ind = c,  palette = "jco" ,label = "ind", geom.var = "point", repel = F, alpha.var ="contrib")
# clustering

cc <- color_Crucero[match(cruise, names(color_Crucero))]
names(cc) <- rownames(ftable)


library(dendextend)
library(gplots)


ftable %>%
  dist(method="euclidean") %>%
  hclust(method="complete") %>% 
  as.phylo.dendrogram() %>%
  plot(type = "phylogram", tip.color = cc, cex = 0.7)

ftable %>%
  dist(method="euclidean") %>%
  hclust(method="complete") %>% 
  as.dendrogram() -> dend


#labels_colors(dend) <- cc
labels(dend) <- s # paste(samples, cruise, sep = '_')


#


pick_top <- function(x, y, top = 10) {

  # x <- vector of abundance
  # y <- vector of taxonomic name
  
  ordered <- order(x, decreasing = T)
  topPos <- head(ordered, top)
  taxPos <- y[topPos]
  
  return(taxPos)
}

v9_fam <- tbl_agglom(tax1, ab1, Level = Level, fkingdom = 'Animalia')

v9_fam %>%
  as_tibble() %>%
  select_at(vars(-contains('X017_X07'))) -> v9_fam

sam <- names(v9_fam)[-1]



apply(v9_fam[-1], 2, pick_top, top = 10, y = v9_fam$Family) %>%
  as_tibble() %>%
  pivot_longer(all_of(sam), names_to = 'id_sample', 
               values_to = Level) %>%
  distinct_at(all_of(Level)) %>%
  inner_join(v9_fam) -> fam_top

# tags
sam <- names(fam_top)[-1]
c <- sapply(strsplit(sam, "[_]"), `[`, 2)
s <- sapply(strsplit(sam, "[_]"), `[`, 3)
t <- substr(s, 1,1)
cc <- paste(c,s, sep = '_')

table(c)

# redoing heatmaps
m <- as(fam_top[-1], 'matrix')
colnames(m) <- cc
rownames(m) <- fam_top$Family


# heatm


total <- sum(rowSums(fam_top[-1]))
tax_sum <- round((rowSums(fam_top[-1]) / total) * 100, 4)


set.seed(20200511)
path1 <- '/Users/cigom/metagenomics/MG_18S/multirun_xiximis/downstream_X4X5X6X7/'

png(paste0(path1,"superheat.png"), 
    height = 11, width = 17, res = 300,
    units = 'in')

apply(m, 2, function(x) x/sum(x) * 100) %>%
  superheat(
          scale = F,
          row.dendrogram = T,
          clustering.method = 'hierarchical',
          dist.method = 'euclidean',
          membership.cols = c,
          left.label.text.size = 4.0,
          left.label.col = 'white',
          # bottom.label.text.size = 2.0,
          # bottom.label.text.angle = 75,
          grid.hline = T,
          grid.vline = T,
          yt = colSums(fam_top[-1]),
          yt.plot.type = "bar",
          yt.axis.name = "Total Reads",
          legend.width = 4,
          legend.text.size = 12,
          legend.breaks = c(0,20,40,60,80))

dev.off()

# select samples
mset1 <- m[,grep('X04', colnames(m))]

# ordenar taxones con respecto a X04 
# superheat(mset1,
#           scale = F,
#           row.dendrogram = T,
#           col.dendrogram = T,
#           clustering.method = 'hierarchical',
#           dist.method = 'euclidean', 
#           print.plot = F) -> sh

# u ordenar con respecto a todo los datos 
superheat(apply(m, 2, function(x) x/sum(x) * 100),
          scale = F,
          row.dendrogram = T,
          col.dendrogram = T,
          clustering.method = 'hierarchical',
          dist.method = 'euclidean',
          print.plot = F) -> sh

# dist <- dist(t(m), method = "euclidean")
tax_hclust <- sh$order.rows
tax_hclust <- rownames(m)[tax_hclust]
#tax_hclust <- rev(tax_hclust)
# sampleOrder <- colnames(otu_table(phyFamilies))[sampleOrder]

# m %>%
#   dist(method="euclidean") %>%
#   hclust(method="complete") %>% 
#   as.phylo.dendrogram() %>%
#   plot(type = "phylogram", tip.color = cc, cex = 0.7)
  


reorder_cormat <- function(matrix){
  # Use correlation between variables as distance
  dd <- dist((1-matrix)/2, method = "euclidean")
  hc <- hclust(dd, method = "centroid")
  matrix_r <- matrix[hc$order, hc$order]
  return(matrix_r)
}


# tl.col.hclust <- color_Crucero[match(cruise_hclust, names(color_Crucero))]
# names(tl.col.hclust) <- sample_hclust

cleanClass <- c('Insecta', 'Mammalia', 'Arachnida')

tax1 %>%
  as_tibble() %>%
  select(Ranks) %>%
  filter(Kingdom == 'Animalia') %>%
  filter(!Class %in%  cleanClass) %>%
  filter(Family != 'Undetermined_R') %>%
  group_by(Family) %>%
  distinct_at(all_of(Level), .keep_all = T) %>%
  select(Phylum,Class,Order,Family) -> tax

plotHeat <- function(mset, facet = TRUE, col_hclust = FALSE) {
  
  if(col_hclust) {
    round(cor(mset), 2) %>% reorder_cormat() %>% rownames() -> sample_hclust
    sample_hclust <- sapply(strsplit(sample_hclust, "[_]"), `[`, 2)
  } else {
    sample_hclust <- colnames(mset)
    sample_hclust <- sapply(strsplit(sample_hclust, "[_]"), `[`, 2)
  }


 p <-  mset %>% 
    as_tibble(rownames = 'Family') %>%
    mutate_at(vars(!all_of(Level)), function(x) {x / sum(x) * 100 }) %>%
    pivot_longer(cols = !all_of(Level), values_to = "ra", 
                 names_to = "station") %>%
    inner_join(tax) %>%
    mutate(Estacion = sapply(strsplit(station, "[_]"), `[`, 2),
           Crucero = sapply(strsplit(station, "[_]"), `[`, 1)) %>%
    mutate(Estacion = factor(Estacion , levels = sample_hclust),
           Family = factor(Family, levels = tax_hclust),
           ra = ifelse(ra == 0, NA, ra)) %>%
    ggplot(aes(x = Estacion, y = Family, fill = ra, color = Crucero))+
    geom_tile(color = "black") +
    ggsci::scale_fill_material("indigo",  
                               name="Abundancia\nRelativa", na.value = 'white',
                               labels = scales::percent_format(scale = 1)) +
    labs(x = NULL, y = NULL) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1),
          axis.text.y = element_text(face = 'bold.italic', size = 12)) +
    guides(fill = guide_colorbar(barheight = unit(9, "in"),
                                 ticks.colour = "black", 
                                 frame.colour = "black",
                                 label.theme = element_text(size = 16))
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, vjust = 1),
      strip.text.y = element_text(
        angle = 0, 
        size = 12),
      strip.background = element_rect(colour = "black", 
                                      fill = "transparent",
                                      size = 0.4),
      panel.spacing = unit(0.007, "lines"))
  
  if(facet) {
    p + facet_grid(Phylum + Class ~ ., 
               scales = "free", space = "free" )
  } else
    return(p)

}

p1 <- plotHeat(mset1)

# 2
mset2 <- m[,grep('X05', colnames(m))]
p2<- plotHeat(mset2)

# 3

mset3 <- m[,grep('X06', colnames(m))]
p3<- plotHeat(mset3)

ggsave(p1, filename = 'heatmap_X4.png', dpi = 300, path = path1,
       units = 'in', width = 17, height = 11)
ggsave(p2, filename = 'heatmap_X5.png', dpi = 300, path = path1,
       units = 'in', width = 17, height = 11)
ggsave(p3, filename = 'heatmap_X6.png', dpi = 300, path = path1,
       units = 'in', width = 17, height = 11)

p1nf <- plotHeat(mset1, facet = F, col_hclust = F) + 
  viridis::scale_fill_viridis()
p2nf <- plotHeat(mset2, facet = F, col_hclust = F) + 
  viridis::scale_fill_viridis()
p3nf <- plotHeat(mset3, facet = F, col_hclust = F) + 
  viridis::scale_fill_viridis()



ggsave(p1nf, filename = 'heatmapnf_X4.png', dpi = 300, path = path1,
       units = 'in', width = 17, height = 11)
ggsave(p2nf, filename = 'heatmapnf_X5.png', dpi = 300, path = path1,
       units = 'in', width = 17, height = 11)
ggsave(p3nf, filename = 'heatmapnf_X6.png', dpi = 300, path = path1,
       units = 'in', width = 17, height = 11)

library(UpSetR)
m %>% 
  as_tibble(rownames = 'Family') %>%
  mutate_at(vars(!all_of(Level)), function(x) {x / sum(x) * 100 }) %>%
  mutate_at(vars(!all_of(Level)), function(x) {ifelse(x > 0, 1, 0) }) %>%
  select(-Family) %>% data.frame(row.names = rownames(m)) %>% #t() %>%
  data.frame() -> binary_m

# Error: vector memory exhausted (limit reached?)
upset(binary_m,
  nsets =  12, #nrow(m),
  #nintersects = NA,
  # Display them from the most numerous intersection to the least
  order.by = "freq",
  line.size = 0.7,
  point.size = 1.7,
  text.scale = 1.3,
  mb.ratio = c(0.30, 0.70),
  empty.intersections = "off")

m %>% 
  as_tibble(rownames = 'Family') %>%
  mutate_at(vars(!all_of(Level)), function(x) {x / sum(x) * 100 }) %>%
  mutate_at(vars(!all_of(Level)), function(x) {ifelse(x > 0, 1, 0) }) %>%
  pivot_longer(cols = !all_of(Level), values_to = "pa", names_to = "id") %>%
  mutate(station = sapply(strsplit(id, "[_]"), `[`, 2),
         cruise = sapply(strsplit(id, "[_]"), `[`, 1)) %>%
  filter(pa > 0) %>%
  group_by(Family, cruise) %>%
  summarise(Freq = sum(pa)) %>%
  ungroup() %>%
  mutate(frac = ifelse(cruise == 'X04', Freq/47,
                       ifelse(cruise == 'X05', Freq/33,
                       ifelse(cruise == 'X06', Freq/42, Freq)))) %>%
  mutate(frac = round(frac, digits = 3)*100 ) %>%
  mutate(Family = factor(Family, levels = tax_hclust)) -> FreqSamples

freqplot <- FreqSamples %>%
  ggplot() +
  geom_tile(aes(x = cruise, y = Family, fill = frac),
    colour = "black",
    size = .3, show.legend = F) +
  geom_text(aes(y = Family, x = cruise, label=as.factor(Freq)), size = 7) +
  ggsci::scale_fill_material("indigo",  
                           name="Presencia\nmuestras", na.value = 'grey',
                           labels = scales::percent_format(scale = 1)) +
  labs(x = NULL, y = NULL) +
  theme_bw() + 
  theme(text = element_text(size = 12),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_line(size = .2),
        panel.border = element_blank(),
        axis.text.y = element_text(face = 'bold.italic'))
  # guides(fill = guide_colorbar(barheight = unit(9, "in"),
  #                              ticks.colour = "black", 
  #                              frame.colour = "black",
  #                              label.theme = element_text(size = 12)))

heatmap <- plotHeat(mset, facet = F, col_hclust = F) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

library(patchwork)
plot <- freqplot + heatmap + plot_layout(widths = c(1, 2))
ggsave(plot, filename = 'patchwork.png', 
       dpi = 300, path = path1,
       units = 'in', width = 5, height = 10)
filename = 'heatmap.png', dpi = 300, path = path1,
units = 'in', width = 17, height = 11)