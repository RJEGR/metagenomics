
library("FactoMineR")
library('FactoInvestigate')
library("factoextra")
library("corrplot")
library(tidyverse)

color_Crucero <- c('X04'='#00AFBB','X05'= '#E7B800', 'X06'='#FC4E07', 'X07'='red')

cleanTax <- c('Arachnida', 'Mammalia', 'Insecta', 'Undetermined_R')

Level <- 'Family'
# agglomerate families
tbl_agglom <- function(tax_tbl, count_tbl, Level = names(tax_tbl)[4]) {
  
  scheck <- identical(rownames(tax_tbl), rownames(count_tbl))
  
  if(!scheck) {
    sortt <- rownames(tax_tbl)
    count_tbl <- count_tbl[match(rownames(count_tbl), sortt),]
  }
  
  tbl <- cbind(tax_tbl, count_tbl)
  sam <- colnames(count_tbl)
  
  tbl_agglom <- tbl %>%
    select(Level, sam) %>%
    drop_na(Level) %>% # to clean the redundant names (NA assignation)
    group_by_at(vars(one_of(Level))) %>%
    summarise_at(vars(sam), sum) %>%
    ungroup()
  
  return(tbl_agglom)
}

 
v9_fam <- tbl_agglom(tax1, ab1, Level = Level)
coi_fam <- tbl_agglom(tax2, ab2, Level = Level)


v9_fam <- v9_fam[v9_fam$Family != 'Undetermined_R',]

total <- sum(rowSums(dat))
tax_sum <- round((rowSums(dat) / total) * 100, 4)

coi_fam %>%
  select(-Level) %>%
  mutate(total = sum(rowSums(.))) %>% 
  mutate(pct = round((rowSums(.) / total) * 100, 4)) %>%
  select(pct)

# transform-data-count, how?
library(ggforce)
sam <- colnames(coi_fam)[-1]
coi_fam %>% 
  filter(!all_of(Level) %in% cleanTax) %>%
  # mutate_at(vars(!all_of(Level)), function(x) {x / sum(x) * 100 }) %>%
  pivot_longer(cols = !all_of(Level), names_to = 'id_sample', 
               values_to = 'ab') %>%
  filter(ab > 1) %>% # Not!!
  mutate(cruise = sapply(strsplit(id_sample, "[.]"), `[`, 1),
         sample = sapply(strsplit(id_sample, "[.]"), `[`, 2),
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

# 18S
# c <- 
# s <- sapply(strsplit(sam, "[_]"), `[`, 3)
# v9_fam
v9_fam %>% 
  filter(!all_of(Level) %in% cleanTax) %>%
  # mutate_at(vars(!all_of(Level)), function(x) {x / sum(x) * 100 }) %>%
  pivot_longer(cols = !all_of(Level), names_to = 'id_sample', 
               values_to = 'ab') %>%
  #filter(ab > 1) %>% # Not!!
  mutate(cruise = sapply(strsplit(id_sample, "[_]"), `[`, 2),
         sample = sapply(strsplit(id_sample, "[_]"), `[`, 3),
         Family = reorder(Family, ab)) -> dv1

coi_fam %>% 
  filter(!all_of(Level) %in% cleanTax) %>%
  # mutate_at(vars(!all_of(Level)), function(x) {x / sum(x) * 100 }) %>%
  pivot_longer(cols = !all_of(Level), names_to = 'id_sample', 
               values_to = 'ab') %>%
  #filter(ab > 1) %>% # Not!!
  mutate(cruise = sapply(strsplit(id_sample, "[.]"), `[`, 1),
         sample = sapply(strsplit(id_sample, "[.]"), `[`, 2),
         Family = reorder(Family, ab)) -> dv2

dv2 %>%
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
# top

v9_fam %>% mutate_at(vars(!all_of(Level)),
                      function(x) {x / sum(x) * 100 }) # -> v9_fam

# select(top)
dv2 %>%
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

# For COI
# c <- sapply(strsplit(sam, "[.]"), `[`, 1)
# s <- sapply(strsplit(sam, "[.]"), `[`, 2)
# t <- substr(s, 1,1)

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


pick_top <- function(x,pos, top = 10) {
  Level <- 'Family'
  
  x <- data.frame(x)
  ordered <- order(x[,pos], decreasing = T)
  topPos <- head(ordered, top)
  taxPos <- taxname[topPos]
  
  return(taxPos)
}

sam <- names(v9_fam)[-1]

apply(v9_fam, 2, pick_top, top = 10) %>%
  as_tibble() %>%
  select(-Level) %>%
  pivot_longer(all_of(sam), names_to = 'id_sample', 
               values_to = Level) %>%
  filter_at(vars(id_sample), all_vars(!grepl('X017_X07',.))) %>%
  #group_by(Level) %>%
  distinct_at(all_of(Level)) %>%
  inner_join(v9_fam) %>%
  select_at(vars(-contains("X017_X07"))) -> fam_top

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

set.seed(20200511)

superheat(m,
          scale = F,
          row.dendrogram = T,
          clustering.method = 'hierarchical',
          dist.method = 'euclidean',
          membership.cols = c,
          left.label.text.size = 2.0,
          left.label.col = 'white',
          # bottom.label.text.size = 2.0,
          # bottom.label.text.angle = 75,
          grid.hline = F,
          grid.vline = T)


# dist <- dist(t(m), method = "euclidean")
# sampleOrder <- 
# sampleOrder <- colnames(otu_table(phyFamilies))[sampleOrder]

# m %>%
#   dist(method="euclidean") %>%
#   hclust(method="complete") %>% 
#   as.phylo.dendrogram() %>%
#   plot(type = "phylogram", tip.color = cc, cex = 0.7)
  
# select samples
mset <- m[,grep('X06', colnames(m))]

reorder_cormat <- function(matrix){
  # Use correlation between variables as distance
  dd <- dist((1-matrix)/2, method = "euclidean")
  hc <- hclust(dd, method = "centroid")
  matrix_r <- matrix[hc$order, hc$order]
  return(matrix_r)
}

round(cor(mset), 2) %>%
  reorder_cormat() %>% rownames() -> sample_hclust

round(cor(t(mset)), 2) %>%
  reorder_cormat() %>% rownames() -> tax_hclust

cruise_hclust <- sapply(strsplit(sample_hclust, "[_]"), `[`, 1)

tl.col.hclust <- color_Crucero[match(cruise_hclust, names(color_Crucero))]
names(tl.col.hclust) <- sample_hclust

# heatmap(as.matrix(dist(t(m))))$colInd
# heatmap(as.matrix(dist(m)))$colInd
tax1 %>%
  select(Ranks) %>%
  as_tibble() %>%
  distinct_at(all_of(Level), .keep_all = T) %>%
  select(Phylum,Class,Order,Family) -> tax
  
m %>% 
  as_tibble(rownames = 'Family') %>%
  mutate_at(vars(!all_of(Level)), function(x) {x / sum(x) * 100 }) %>%
  pivot_longer(cols = all_of(sample_hclust)) %>%
  mutate(Estacion = sapply(strsplit(name, "[_]"), `[`, 2),
         Crucero = sapply(strsplit(name, "[_]"), `[`, 1)) %>%
  mutate(name = factor(name , levels = sample_hclust),
         Family = factor(Family, levels = tax_hclust),
         value = ifelse(value == 0, NA, value)) %>%
  inner_join(tax) %>%
  ggplot(aes(x = Estacion, y = reorder(Family, value), fill = value, color = Crucero))+
  geom_tile(color = "grey")+
  ggsci::scale_fill_material("indigo",  name="Relative\nAbundance") +
  labs(x = "", y = "") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1),
        axis.text.y = element_text(face = 'bold.italic')) +
  guides(fill = guide_colorbar(barheight = unit(7, "cm"),  
                               ticks.colour = "black", frame.colour = "black")) +
    facet_grid(Phylum+Class ~ ., 
               scales = "free", space = "free" ) + 
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
    )

