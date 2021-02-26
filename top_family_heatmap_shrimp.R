
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

obj %>%
  select_if(is.double) %>%
  names() -> colNames


# clean taxonomy

obj %>%
  mutate_at(ranks, funs(str_replace_all(., c("_[1-9]" = "", ".+_unclassified"=NA_character_, "_"=" ", "Unclassified"=NA_character_, "NA"=NA_character_, "sp."=NA_character_, "Undetermined"=NA_character_, "Unknown"=NA_character_)))) -> obj

agglom_lev <- "Family"


obj %>% distinct_at(agglom_lev)

# prepare data to select top

obj %>%
  select_at(vars(c(agglom_lev, colNames))) %>%
  rename("Level" = agglom_lev) %>%
  # filter(grepl('ceae', Level)) %>% # keep all 
  group_by(Level) %>%
  summarise_at(vars(colNames), sum) %>%
  ungroup() -> agg_wide

# test prevalence/ab

prevelancedf = apply(X = agg_wide[-1],
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

df <- data.frame(Prevalence = prevelancedf, 
                TotalAbundance = rowSums(agg_wide[-1]),
                Family = agg_wide$Level) %>% 
  left_join(obj %>% select_at(ranks[1:5]) %>% distinct_at(agglom_lev, .keep_all = TRUE)) %>%
  as_tibble() %>%
  rename("Level" = agglom_lev) %>%
  arrange(desc(TotalAbundance))


pd <- position_dodge(0.1)

df %>% group_by(Phylum) %>% summarise(n = sum(TotalAbundance), 
                                      mean = mean(Prevalence),
                                      N = sum(!is.na(Phylum)),
                                      sd = sd(Prevalence, na.rm = T),
                                      se = sd  / sqrt(N)) %>% 
  arrange(desc(n)) -> summary_prevalence

summary_prevalence %>%
  ggplot(aes(n, mean)) +
  # geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1, position = pd) +
  # geom_line(position = pd) +
  # geom_point(aes(size = N), alpha = 0.7, position = pd) +
  scale_x_log10() +
  ggrepel::geom_text_repel(aes(label = Phylum, color = N), size = 4) +
  labs(y = "Prevalence (mean)", x = "Total Abundance")

keepPhyla <- summary_prevalence$Phylum[summary_prevalence$n/sum(summary_prevalence$n) >= 0.001]

df %>%
  # filter(Phylum %in% keepPhyla) %>%
  mutate(Level = ifelse(is.na(Level), "Incomplete", "Complete")) %>%
  mutate(Phylum= ifelse(Phylum %in% keepPhyla, Phylum, "Low Taxa")) %>%
  mutate(Phylum = factor(Phylum, levels = c(keepPhyla, "Low Taxa"))) %>%
  # mutate(Top = ifelse(Level %in% fam_top$Level, TRUE, FALSE)) %>%
  mutate(Prevalence = Prevalence/length(colNames)) %>%
  ggplot(aes(TotalAbundance, Prevalence)) + 
  geom_point(size = 2, alpha = 0.7, aes(color = Level)) + 
  geom_hline(yintercept = 0.10, alpha = 0.5, linetype = 2) +
  scale_x_log10() +
  labs(y = "Prevalence (Frac. Samples)", x = "Total Abundance (log10)", color = "Lineage") +
  facet_wrap(~Phylum) +
  theme_bw(base_size = 17) +
  theme(
    legend.position = "top",
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12)) -> saveP

ggsave(saveP, filename = "Family_prevalence.png", path = dir, 
       width = 8, height = 8)

save(df, keepPhyla, file = paste0(dir, "ANCOM_BC_inputs.RData"))

pick_top <- function(x, y, top = 10) {
  
  # x <- vector of abundance
  # y <- vector of taxonomic name
  
  ordered <- order(x, decreasing = T)
  topPos <- head(ordered, top)
  taxPos <- y[topPos]
  
  return(taxPos)
}

obj %>%
  select_at(vars(ranks)) %>%
  distinct_at(agglom_lev, .keep_all = T) %>%
  rename("Level" = agglom_lev) %>%
  left_join(agg_wide) %>%
  filter(Phylum %in% keepPhyla) %>%
  select_at(vars(c("Level",colNames))) %>%
  mutate(Level = ifelse(is.na(Level), "Incomplete", Level)) -> agg_wide_filtered


apply(agg_wide_filtered[-1], 2, pick_top, top = 10, 
      y = agg_wide_filtered$Level) %>%
  as_tibble() %>%
  pivot_longer(all_of(colNames), names_to = 'Index', 
               values_to = "Level") %>%
  distinct(Level) %>%
  inner_join(agg_wide_filtered) -> fam_top




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
  mutate(Family = ifelse(is.na(Family), "Incomplete", Family)) %>%
  select_at(vars(ranks)) %>%
  distinct_at(agglom_lev, .keep_all = T) %>%
  rename("Level" = agglom_lev) %>%
  inner_join(fam_top, by = "Level") %>%
  mutate_at(colNames, raf) %>%
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ra') %>%
  filter(ra > 0) %>% inner_join(mtd) -> dataHeat

# sanity check ----



# set left-panel (Phylum) ordering 
dataHeat %>% group_by(Phylum) %>% summarise(t = sum(ra)) %>% arrange(desc(t)) %>% pull(Phylum) -> PhylumLevel

TimeLev <- c("Farm", 0,20,40,60,80)

dataHeat %>%
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) %>%
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

save(obj, mtd, colNames, keepPhyla,
     tax_hclust, ranks, file = paste0(dir, "objects.RData"))

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


