# Set colors
color_Crucero <- c('X04'='#8c510a','X05'= '#bf812d', 'X06'='#dfc27d', 'X07'='red')

color_strategy <- c('OTU'='#253494', 'ASV'='#636363')

Ranks <- c("Kingdom", "Phylum", "Class", "Order", 
           "Family", "Genus", "Species")

# Set paths ----
p1 <- '/Users/cigom/metagenomics/MG_18S/multirun_xiximis/downstream_X4X5X6X7/'

p2 <- '/Users/cigom/metagenomics/MG_COI/multirun_20200507_COI/'

out_path <- '~/metagenomics/'

# Count data ----

c1 <- dir(path = p1, 
          pattern = '*list.shared', 
          full.names = T)

ab <- read.table(c1)

# metadata ----
m1 <- data.frame(
  row.names = names(ab),
  Crucero = sapply(strsplit(names(ab), "[_]"), `[`, 2),
  Station =  sapply(strsplit(names(ab), "[_]"), `[`, 3),
  Strategy = '18S')

# taxonomy ----
t1 <- dir(path = p1, 
          pattern = '*cons.taxonomy', 
          full.names = T)

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

dim(tax <- read_rdp2(t1, header = T))

# Make phyloseq ----

library(phyloseq)

identical(names(ab),rownames(m1))
identical(rownames(ab), rownames(tax))

physeq <- phyloseq(otu_table(ab, taxa_are_rows = TRUE),
                   tax_table(as(tax, 'matrix')), 
                    sample_data(m1))

# Removing any abundance of zero and

#cleanClass <- c('Arachnida', 'Mammalia', 'Insecta','Reptilia', 'Aves')

cleanTax <- 'Undetermined_R'


physeq %>% 
  subset_taxa(Kingdom == "Animalia") %>%
  subset_taxa(!Family %in%  cleanTax) %>%
  prune_taxa(taxa_sums(.) > 0,.) -> physeq

# count the number of families
taxrank <- 'Family'
top <- 50
family <- tax_glom(physeq, taxrank = taxrank, NArm = FALSE)

otu_table <- as(otu_table(family), 'matrix') %>% as_tibble()
tax_table <- as(tax_table(family), 'matrix') %>% as_tibble() %>% select(all_of(taxrank))

dat <- data.frame(otu_table, row.names = tax_table$Family)

total <- sum(rowSums(dat))
tax_sum <- round((rowSums(dat) / total) * 100, 4)

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

# plot(ecdf(yr_plot_clus))
# plot(ecdf(yr_plot_max))

sam <- names(dat)

dat %>% 
  filter(rownames(.) %in% names(yr_plot_max)) %>%
  mutate_at(vars(sam),
          function(x) {x / sum(x) * 100 }) -> select_top_pct 
  
rownames(select_top_pct) <- names(yr_plot_max)
superheat(select_top_pct,
          col.dendrogram = T,
          clustering.method = 'hierarchical',
          linkage.method = 'complete',
          dist.method = 'euclidean',
          pretty.order.rows = T,
          yr = yr_plot_max,
          yr.plot.type = 'scatterline',
          left.label.text.size = 2.0,
          left.label.col = 'white',
          bottom.label = 'variable',
          bottom.label.text.size = 2.0,
          bottom.label.text.angle = 75,
          grid.hline = F,
          grid.vline = F)

# choose top of families
topFamilies <- names(yr_plot_max)

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

# mdf <- agglom_rank(physeq, taxrank, fra)
mdf <- phyFamilies %>%
  transform_sample_counts(function(x) 1E2 * x/sum(x)) %>%
  psmelt()


mdf[mdf$Abundance == 0, 'Abundance'] <- NA


colourCount <- length(unique(mdf$Class)) 
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
fvalues = getPalette(colourCount)

plot1 <- mdf %>% 
  filter(Crucero %in% c('X04', 'X06') ) %>%
  ggplot(aes(x=Station, 
                   y=Family, 
                   color=Class)) +
  geom_point(aes(size=Abundance))+ facet_grid(.~Crucero, scales="free_x") +
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


# heatmap

phyFamilies %>%
  transform_sample_counts(function(x) 1E2 * x/sum(x)) -> heatFam

dist <- distance(heatFam, method = "euclidean")
sampleOrder <- heatmap(as.matrix(dist))$colInd
sampleOrder <- colnames(otu_table(heatFam))[sampleOrder]

taxaOrder = names(sort(taxa_sums(heatFam), decreasing = TRUE))



otu_table(heatFam)[otu_table(heatFam) == 0] <- NA

#taxaOrder = rownames(sort(unique(tax_table(Family)[,'Orden']), decreasing=TRUE))
colourCount <- 8 
getPalette = colorRampPalette(brewer.pal(colourCount, "YlGnBu"))
fvalues = getPalette(colourCount)

h <- NULL
heatFamPlot <- subset_samples(heatFam, Crucero %in% c('X05', 'X07'))
h <- plot_heatmap(heatFamPlot, 
                  sample.label = "Station", 
                  sample.order = sampleOrder,
                  taxa.label = taxrank,
                  taxa.order = taxaOrder,
                  na.value = "white", trans = NULL,
                  low = fvalues[1:4],
                  high = fvalues[4:8]
) +
  labs(fill = "Abundancia\nRelativa (%)") +
  theme(legend.position = "left")

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


plot
