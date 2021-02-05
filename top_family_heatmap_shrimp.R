
rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)
mtd_file <- list.files(path = dir, pattern = "CIAD_MappingFile.txt$", full.names = TRUE)

library(readxl)
library(tidyverse)


obj <- read_xlsx(files[2])

mtd <- read.csv(mtd_file, sep = "\t") %>% 
  select(Sample_IDs, Tissue, Time) %>%
  rename(Index = Sample_IDs)

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")


# clean taxonomy

obj %>%
  mutate_at(ranks, 
            funs(str_replace_all(., c("_[1-9]" = "")))) -> obj

obj %>%
  select_if(is.double) %>%
  names() -> colNames

agglom_lev <- "Family"

# prepare data to select top

obj %>%
  select_at(vars(c(agglom_lev, colNames))) %>%
  rename("Level" = agglom_lev) %>%
  # mutate(Level = str_replace_all(Level, c("_[1-9]" = ""))) %>%
  filter(grepl('ceae', Level)) %>%
  group_by(Level) %>%
  summarise_at(vars(colNames), sum) %>%
  ungroup() -> agg_wide


pick_top <- function(x, y, top = 10) {
  
  # x <- vector of abundance
  # y <- vector of taxonomic name
  
  ordered <- order(x, decreasing = T)
  topPos <- head(ordered, top)
  taxPos <- y[topPos]
  
  return(taxPos)
}

apply(agg_wide[-1], 2, pick_top, top = 10, 
      y = agg_wide$Level) %>%
  as_tibble() %>%
  pivot_longer(all_of(colNames), names_to = 'Index', 
               values_to = "Level") %>%
  distinct(Level) %>%
  inner_join(agg_wide) -> fam_top

# make tax clustering 

library(superheat)

m <- data.frame(fam_top[-1])
rownames(m) <- fam_top$Level

raf <- function(x) x/sum(x) * 100

superheat(apply(m, 2, raf),
          scale = F,
          row.dendrogram = T,
          col.dendrogram = T,
          clustering.method = 'hierarchical',
          dist.method = 'euclidean',
          print.plot = F) -> sh

tax_hclust <- sh$order.rows
tax_hclust <- rownames(m)[tax_hclust]

obj %>%
  select_at(vars(ranks)) %>%
  distinct_at(agglom_lev, .keep_all = T) %>%
  rename("Level" = agglom_lev) %>%
  inner_join(fam_top, by = "Level") %>%
  mutate_at(colNames, raf) %>%
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ra') %>%
  filter(ra > 0) %>% inner_join(mtd) -> dataHeat

# sanity check ----

dataHeat %>% group_by(Tissue, Index) %>% summarise(sum(ra))

# set left-panel (Phylum) ordering 
dataHeat %>% group_by(Phylum) %>% summarise(t = sum(ra)) %>% arrange(desc(t)) %>% pull(Phylum) -> PhylumLevel

dataHeat %>%
  mutate(Level = factor(Level, levels = tax_hclust),
         Phylum = factor(Phylum, levels = PhylumLevel),
         ra = ifelse(ra == 0, NA, ra)) %>%
  ggplot(aes(y = Level, x = Index, fill = ra)) +
  geom_tile() +
  # facet_grid( ~ as.factor(Tissue) , scales = "free", space = "free", switch = "x") +
  geom_tile(color = "black") +
  ggsci::scale_fill_material("blue-grey",  
                             name="Relative\nAbundance\n", 
                             na.value = 'white',
                             limits = c(0,100),
                             labels = scales::percent_format(scale = 1)) +
  # ggh4x::facet_nested(~ Tissue + Time, scales = "free", space = "free") +
  ggh4x::facet_nested(Phylum + Class ~ Tissue + Time,
                      scales = "free", space = "free" ) +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_family = "GillSans", base_size = 14) +
  guides(fill = guide_colorbar(barheight = unit(9, "in"), 
                               barwidth = unit(0.5, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 14))
  ) +
  
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text.y = element_text(
      angle = 0, 
      size = 14),
    strip.background = element_rect(colour = "black", 
                                    fill = "transparent",
                                    size = 0.4),
    panel.spacing = unit(0.007, "lines")) -> heatPlot

ggsave(heatPlot, filename = "heatmap.png", path = dir, 
       width = 22, height = 16)

save(obj, mtd, colNames, tax_hclust, ranks, file = paste0(dir, "objects.RData"))

# test features diversity ----

source("~/Documents/GitHub/metagenomics/estimate_richness.R")

obj %>%
  select_at(vars(c(colNames))) %>%
  estimate_richness(., measures = c("Observed", "Shannon", "Chao1")) %>%
  as_tibble(rownames = "Index") %>%
  inner_join(mtd) -> diversityDat

diversityDat %>%
  ggplot(aes(x = Chao1, y = Observed)) +
  geom_point(data = diversityDat %>% select(-Time), 
             colour = "grey70", alpha = 0.3, size = 2.5) +
  geom_point(aes(colour = Tissue), size = 2.5) +
  xlab("Predicted (Chao1 Estimator)") +
  ylab("Observed") +
  facet_wrap(~ Time ) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  ggsci::scale_color_rickandmorty() -> div1 

ggsave(div1, filename = "obs_vs_chao_tissue_time.png", path = dir, 
       width = 12, height = 6)

diversityDat %>%
  ggplot(aes(x = Chao1, y = Observed)) +
  geom_point(data = diversityDat %>% select(-Tissue), colour = "grey70",
             alpha = 0.3) +
  geom_point(aes(colour = Shannon)) +
  facet_wrap(~ Tissue ) +
  ggrepel::geom_text_repel(max.overlaps = 30, size = 2.5,
                           aes(Chao1, Observed,
                               label = Index, color = Shannon)) +
  labs(xlab = "Predicted (Chao1 Estimator)",
       ylab = "Observed") +
  scale_colour_viridis_c(option = "magma", begin = 0.2, end = 0.8) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "top") +
  guides(color = guide_colorbar(barheight = unit(0.2, "in"),
                                barwidth = unit(8, "in"),ticks.colour = "black", frame.colour = "black",label.theme = element_text(size = 10))) -> div2

ggsave(div2, filename = "obs_vs_chao_tissue_shannon.png", path = dir, 
       width = 12, height = 6)

diversityDat %>%
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  # mutate(Time = ifelse(Time %in% "Farm", 0, Time)) %>%
  # mutate(Time = as.double(Time)) %>%
  # mutate(Time = factor(Time, levels = c(0,20,40,60,80))) %>%
  ggplot(aes(x = Time, y = Shannon)) +
  geom_boxplot(lwd = 0.5) +
  # stat_boxplot(geom ='errorbar', linetype = "dotted" ) +
  geom_point(aes(x = Time, y = Shannon, color = Tissue)) +
  theme_classic(base_size = 16, base_family = "GillSans") +
  ggsci::scale_color_rickandmorty() -> div3
div3

ggsave(div3, filename = "Time_shannon_tissue.png", path = dir, 
       width = 10, height = 4)

# betadiv ----

library(vegan)

# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html

beta_mds <- function(df) {
  df %>%
    t() %>%
    vegan::vegdist(., method="bray") %>%
    vegan::metaMDS(.) -> mds
  
  mds_data <- as.data.frame(mds$points)
  mds_data$Index <- rownames(mds_data)
  
  return(mds_data)
  
}

obj %>%
  select_at(vars(c(colNames))) %>%
  beta_mds(.) %>%
  left_join(., mtd) %>%
  ggplot(., aes(x = MDS1, y = MDS2, fill = Tissue, color = Tissue)) +
  ggforce::geom_mark_ellipse(aes(label = Tissue, group = Tissue)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(x = MDS1, y = MDS2, label = Index, 
                                max.overlaps = 50), color = "black", fill = "transparent") +
  ggsci::scale_color_rickandmorty() +
  ggsci::scale_fill_rickandmorty() +
  labs(caption = "Multidimensional Scaling (NMDS) of\nPairwise beta diversity indexes (Bray-Curtis index)\n ") +
  theme_light(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "none") -> Pbeta

ggsave(Pbeta, filename = "Bray_curtis_NMDS.png", path = dir, 
       width = 12, height = 10)


# metacoder ----

# install.packages("devtools")
# devtools::install_github("grunwaldlab/metacoder")

library(metacoder)

obj %>%
  select_at(vars(ranks)) %>%
  distinct(Family, .keep_all = T) %>%
  inner_join(fam_top, by = Level) %>%
  parse_tax_data(class_cols = ranks, named_by_rank = TRUE) -> x

# getting per-taxon information
x$data$tax_abund <- calc_taxon_abund(x, "tax_data",
                                     cols = mtd$Index)

x$data$diff_table <- compare_groups(x, 
                                    dataset = "tax_abund",
                                    # What columns of sample data to use
                                    cols = mtd$Index,
                                    # What category each sample is assigned to
                                    groups = mtd$Tissue) 

set.seed(1)

heat_tree_matrix(x,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "differential_heat_tree_Shrimp.pdf")

# or


x$data$tax_occ <- calc_n_samples(x, "tax_abund", groups = mtd$Tissue, cols = mtd$Index)

set.seed(1) 
heat_tree(x, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Intestine, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford",
          output_file = "Intestine_heat_tree_Shrimp.pdf") 

heat_tree(x, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Hepatopancreas, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford",
          output_file = "Hepatopancreas_heat_tree_Shrimp.pdf") 

heat_tree(x, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Stomach, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford",
          output_file = "Stomach_heat_tree_Shrimp.pdf") 

# pos-clustering w/ Lulu ----

# devtools::install_github("tobiasgf/lulu") 


vsearch.file <- list.files(path = dir, pattern = "match_list.txt", full.names = T)
matchlist <- read.table(vsearch.file , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

feature_tab <- obj %>% select_if(is.double) %>% 
  data.frame(., row.names = obj$`Feature ID`)

# vegan::rarecurve(t(feature_tab))

# # # # # # # # # # #
# Run lulu algorithm
# # # # # # # # # # # 

library(lulu)

min_r <- 1
min_match <- 98
min_cooccurence <- 0.95 

curated_result <- lulu(feature_tab, matchlist, minimum_ratio_type = "min", 
                       minimum_ratio = min_r, 
                       minimum_match = min_match, 
                       minimum_relative_cooccurence = min_cooccurence)


otu_map <- curated_result$otu_map
parent <- otu_map[otu_map$curated == "parent",]

cat("\n...Number of Processed Features are:", length(otu_map$curated), "\n")
cat("Curated in ~real biological Features (parents):",table(otu_map$curated)[2], "\n")
cat("And", table(otu_map$curated)[1], "daughters (merged) Features\n")

curated_table <- curated_result$curated_table

tax <- obj %>% select_at(ranks) %>%
  data.frame(row.names = obj$`Feature ID`) %>%
  filter(row.names(.) %in% rownames(curated_table))

testing <- identical(rownames(tax), rownames(curated_table))
cat("\n.... Taxonomy and count table has been curated:",testing, "\n")

# test beta div
curated_table %>%
  beta_mds(.) -> betaLulu


obj %>% # non_curated matrix
  select_at(vars(c(colNames))) %>%
  beta_mds(.) %>%
  left_join(., mtd) %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  geom_point(data = betaLulu, colour = "grey70", alpha = 0.3, size = 2.5) +
  geom_point(aes(fill = Tissue),shape = 21, colour = "grey33") +
  # ggforce::geom_mark_ellipse(aes(label = Tissue, group = Tissue)) +
  # ggsci::scale_color_rickandmorty() +
  ggsci::scale_fill_rickandmorty() +
  labs(caption = "Multidimensional Scaling (NMDS) of\nPairwise beta diversity indexes (Bray-Curtis index)\n ") +
  theme_light(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "top") -> Pbeta

ggsave(Pbeta, filename = "Bray_curtis_NMDS_lulu.png", path = dir, 
       width = 7, height = 5)

# and alpha for lulu
TimeLev <- c("Farm", 0,20,40,60,80)

diversityDat %>% 
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) -> diversityDat

curated_table %>%
  estimate_richness(., measures = c("Observed", "Shannon", "Chao1")) %>%
  as_tibble(rownames = "Index") %>%
  left_join(diversityDat , by = "Index", suffix = c("_Curated", "_Raw")) %>%
  pivot_longer(cols = contains(c("Observed", "Shannon", "Chao1"))) %>%
  separate(col = name, c("richness", "group"), "_") %>%
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) %>%
  filter(richness %in% c("Shannon", "Observed")) -> df_alfa

df_alfa %>% 
  ggplot(aes(x = Tissue, y = value)) +
  facet_grid(richness ~ Time, scales = "free") +
  geom_bar(aes(fill = group), position="dodge", stat="identity") +
  # geom_boxplot(aes(color = group)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1, vjust = 1, size = 10)) +
  ylab("Richness") +
  scale_fill_manual("", values = c("#3182bd", "#de2d26")) -> psave

ggsave(psave, filename = "richness_raw_and_lulu.png", path = dir, 
       width = 7, height = 5)

# are statistically significant different? compare means!!

my_comparisons <- list( c("Intestine", "Stomach") )

library(ggpubr)

psave + 
  stat_compare_means(group.by = "group",paired = TRUE, label = "p.signif", label.y.npc = "top")

source("~/Documents/GitHub/Estadistica_UABC/anova_and_gaussianity.R")

df_alfa %>%
  filter(richness %in% "Shannon") %>%
  # group_by(Tissue, Time, group) %>%
  rename("x" = value, "g" = Index) %>%
  # mutate(g = paste0(Tissue, group, collapse = "_")) %>%
  is_parametric() %>% view()

# por tanto, prueba t para dos muestras dependientes (raw y curated)
# clean NA to unknown


df_alfa %>%
  filter(richness %in% "Shannon") %>%
  group_by(Tissue, group) %>%
  summarise(
    count = n(),
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
  )

df_alfa %>%
  filter(richness %in% "Shannon") %>%
  wilcox.test(value ~ group, data = ., paired = TRUE)

df_alfa %>%
  filter(richness %in% "Shannon") %>%
  pivot_wider(names_from = group) %>%
  mutate(x = Curated-Raw) %>%
  rename("g" = Tissue) -> Dif
  
is_parametric(Dif)
n <- nrow(Dif)
calculate_r_critical_coeff(n) 
(mean(Dif$x))/(sd(Dif$x)/sqrt(n))

qt(p= 0.975, df= n-1, lower.tail=F)


for(i in ncol(tax):1) {
  # it will replace 'unknown_unclassified' terms for 'Unknown' term
  tax[,i] <- sub("unknown", "Unknown", tax[,i])
  tax[,i] <- sub("NA.*", "Unknown", tax[,i])
  #   Sub taxa with *nclassified with Undetermined term
  tax[,i] <- sub(".*nclassified.*", 
                 "Undetermined", tax[,i], perl=TRUE)
  #   Instead of Undetermined classified , use NA
  #tax[,i][grep("Undetermined|Unknown|NA",
  tax[,i][grep("Undetermined|NA",
               tax[,i], perl=TRUE)] <- NA
}


otu_richness <- colSums(curated_table>0)
total_otu <- nrow(curated_table)
betadiversity <- total_otu[i]/mean(otu_richness)

