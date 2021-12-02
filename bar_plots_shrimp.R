rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)
mtd_file <- list.files(path = dir, pattern = "CIAD_MappingFile.txt$", full.names = TRUE)

library(readxl)
library(tidyverse)
library(ggsci)
library(ggsci)

obj <- read_xlsx(files[2])

obj %>%
  select_if(is.double) %>%
  names() -> colNames

mtd <- read.csv(mtd_file, sep = "\t") %>% 
  select(Sample_IDs, Tissue, Time) %>%
  rename(Index = Sample_IDs)

sam <- data.frame(mtd, row.names = mtd$Index) %>% 
  arrange(match(Index, colNames)) %>% mutate_if(is.character, as.factor)

# foregut (stomach), midgut (hepatopancreas), and hindgut (intestine), 

sam %>% mutate(Tissue = recode_factor(Tissue, Intestine = "Hindgut", Hepatopancreas = "Midgut", Stomach = "Foregut")) -> sam

TimeLev <- c("Farm", 0,20,40,60,80)

sam %>% mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) -> mtd

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")


# clean taxonomy

obj %>%
  mutate_at(ranks, funs(str_replace_all(., c("_[1-9]" = "", ".+_unclassified"=NA_character_, "_"=" ", "Unclassified"=NA_character_, "NA"=NA_character_, "sp."=NA_character_, "Undetermined"=NA_character_, "Unknown"=NA_character_)))) -> obj

obj %>%
  mutate_at(ranks, funs(str_replace_na)) %>%
  mutate_at(ranks, funs(str_replace_all(., c("NA" = "")))) -> obj


barTax <- function(obj, colNames, agglom_lev = "Phylum", low_ab = 1) {
  
  obj %>%
    select_at(vars(ranks)) %>%
    rename( "Level" = agglom_lev) %>%
    distinct(Level, .keep_all = T) %>% 
    drop_na(Level) -> tax 
  
  obj %>% 
    pivot_longer(cols = colNames, names_to = "Index", values_to = 'ab') %>%
    filter(ab > 0) %>%
    inner_join(mtd) %>%
    rename( "Level" = agglom_lev) %>%
    group_by(Level, Tissue) %>%
    summarise(ab = sum(ab), Freq = length(Level > 0)) %>%
    group_by(Tissue) %>%
    mutate(RA = (ab / sum(ab)) * 100) %>%
    left_join(tax) %>%
    mutate(Level = ifelse(RA < low_ab, "ZLow", Level)) -> dataViz
  
  labels <- dataViz %>% pull(Level) %>% unique() %>% sort()
  colourCount = length(labels)
  
  if(colourCount > 7) {
    getPalette <- colorRampPalette(pal_uchicago()(7))(colourCount)
  } else
    getPalette <- pal_uchicago(colourCount)
  
  # pal_locuszoom
  # 
  
  getPalette[length(getPalette)] <- "Black"
  labels[length(labels)] <- "Low abundance"
  
  dataViz %>%
    ggplot(aes(x = Tissue, y = RA, fill = Level)) +
    geom_col() +
    coord_flip() +
    labs(x = "Tissue", y = "Relative Abundance (%)", 
         caption = paste0("Low Abundance (Relative Ab <", low_ab, " %)")) +
    scale_fill_manual(agglom_lev, labels = labels, values = getPalette) +
    theme_classic(base_size = 18, base_family = "GillSans")
  
}

PPhylum <- barTax(obj, colNames, "Phylum", low_ab = 1)
PClass <- barTax(obj, colNames, "Class", low_ab = 1)
POrder <- barTax(obj, colNames, "Order", low_ab = 2)

# ph <- unique(PPhylum$data$Level)[-1]
# POrder$data <- POrder$data %>% filter(Phylum %in% ph)
# POrder + facet_wrap(~Phylum, scales = 'free') -> POrder

# ggsave(POrder, filename = "Bar_Order_facet.png", path = dir, 
#        width = 14, height = 7)

ggsave(PPhylum, filename = "Bar_Phylum.png", path = dir, 
       width = 7, height = 3.5)
ggsave(PClass, filename = "Bar_Class.png", path = dir, 
       width = 7, height = 3.5)
ggsave(POrder, filename = "Bar_Order.png", path = dir,
      width = 7.5, height = 4.7)

# by time

facet_bar <- function(obj, colNames, agglom_lev = "Phylum", low_ab = 1) {
  
  TimeLev <- c("Farm", 0,20,40,60,80)
  
  obj %>% 
    pivot_longer(cols = colNames, names_to = "Index", values_to = 'ab') %>%
    filter(ab > 0) %>%
    inner_join(mtd) %>%
    rename( "Level" = agglom_lev) %>%
    group_by(Level, Tissue, Time) %>%
    summarise(Freq = sum(ab > 0), ab = sum(ab)) %>%
    group_by(Tissue, Time) %>%
    mutate(ra = (ab / sum(ab)) * 100) %>%
    mutate(Level = ifelse(ra < low_ab, "ZLow", Level)) %>%
    mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
    mutate(Time = factor(Time, levels = TimeLev)) -> dataViz
  
  dataViz %>% group_by(Time, Tissue) %>% summarise(sum(ra))
  
  labels <- dataViz %>% pull(Level) %>% unique() %>% sort()
  colourCount = length(labels)
  
  if(colourCount > 7) {
    getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
  } else
    getPalette <- pal_locuszoom()(colourCount)
  
  
  getPalette[length(getPalette)] <- "Black"
  labels[length(labels)] <- "Low abundance"
  
  dataViz %>%
    drop_na()%>%
    ggplot(aes(x = Time, y = ra, fill = Level)) +
    geom_col() + # scale_x_reverse() +
    coord_flip() +
    facet_grid(Tissue ~., scales = "free_y") +
    labs(x = "Time", y = "Relative Abundance (%)", 
         caption = paste0("Low Abundance (Relative Ab <", low_ab, " %)")) +
    scale_fill_manual(agglom_lev, labels = labels, values = getPalette) +
    theme_classic(base_size = 18, base_family = "GillSans")
}

PPhylum <- facet_bar(obj, colNames, "Phylum")
PClass <- facet_bar(obj, colNames, "Class", low_ab = 1)
POrder <- facet_bar(obj, colNames, "Order", low_ab = 2)

ggsave(PPhylum, filename = "Bar_Phylum_time.png", path = dir, 
       width = 10, height = 8)
ggsave(PClass, filename = "Bar_Class_time.png", path = dir, 
       width = 10, height = 8)
ggsave(POrder, filename = "Bar_Order_time.png", path = dir, 
       width = 10, height = 8)

# alluvial/metacoder/treemap test ----

obj %>% 
  pivot_longer(cols = colNames, names_to = "Index", values_to = 'ab') %>%
  filter(ab > 0) %>%
  inner_join(mtd) %>%
  # rename( "feature" = `Feature ID`) %>%
  group_by(`Feature ID`, Tissue) %>%
  summarise(ab = sum(ab)) %>%
  pivot_wider(names_from = Tissue, values_from = ab, values_fill = 0) %>%
  left_join(obj %>% select_at(vars(-colNames)), by = "Feature ID") %>%
  ungroup() -> datviz


# metacoder
library(metacoder)

obj  %>% # datviz
  parse_tax_data(class_cols = ranks, named_by_rank = TRUE) -> x


# getting per-taxon information

# cols = c('Hindgut', 'Foregut', 'Midgut')
cols = mtd$Index

x$data$tax_abund <- calc_taxon_abund(x, "tax_data", cols = cols)

x$data$diff_table <- compare_groups(x, 
                                    data = "tax_abund", 
                                    cols = cols, 
                                    groups = mtd$Tissue)


set.seed(1)

output_file <- paste0(dir, "differential_heat_tree_tax_abund.pdf") # png

x <- mutate_obs(x, "diff_table",
                wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"),
                log2_median_ratio = ifelse(wilcox_p_value < 0.05 | is.na(wilcox_p_value), log2_median_ratio, 0))

x %>% 
  mutate_obs("tax_abund", abundance = rowMeans(x$data$tax_abund[cols])) %>%
  filter_taxa(abundance >= 0.001, reassign_obs = c(diff_table = FALSE)) %>%
  heat_tree_matrix(.,
                   data = "diff_table",
                   node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                   node_label = taxon_names,
                   node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                   node_color_range = diverging_palette(), # The built-in palette for diverging data
                   node_color_trans = "linear", # The default is scaled by circle area
                   node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                   edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                   node_size_axis_label = "Number of Features",
                   node_color_axis_label = "Log2 ratio median proportions",
                   layout = "davidson-harel", # The primary layout algorithm
                   initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                   output_file = output_file)


# or by otu-abundance

# Convert counts to proportions
x$data$otu_table <- calc_obs_props(x, data = "tax_data", cols = mtd$Index)

# Get per-taxon counts
x$data$tax_table <- calc_taxon_abund(x, data = "otu_table", cols = mtd$Index)

x$data$diff_table <- compare_groups(x, 
                                    data = "otu_table", 
                                    cols = cols, 
                                    groups = mtd$Tissue)

set.seed(1)

output_file <- paste0(dir, "differential_heat_tree_otu_table.png")

x <- mutate_obs(x, "diff_table",
                wilcox_p_value = p.adjust(wilcox_p_value, method = "fdr"),
                log2_median_ratio = ifelse(wilcox_p_value < 0.05 | is.na(wilcox_p_value), log2_median_ratio, 0))
# not difference w/ otu_ab


# ggalluv (not working )

datviz %>% distinct(Family, .keep_all = T) %>% 
  select_at(vars(ranks[2:6])) %>%
  drop_na(Family) -> tax

datviz %>%
  group_by(Family) %>%
  summarise_at(vars(c('Hindgut', 'Foregut', 'Midgut')), sum) %>%
  drop_na() %>%
  left_join(tax) -> datviz

# library(treemap)
# library(RColorBrewer)
library(ggalluvial)

datviz %>%
  pivot_longer(cols = c('Hindgut', 'Foregut', 'Midgut'), 
               names_to = 'sam', values_to = 'ab') %>%
  filter(ab > 0) %>%
  to_lodes_form(datviz, key = 'sam',
                axes = ranks[2:6]) -> longer_df

longer_df %>% 
  ggplot(aes(x = sam, stratum = stratum, alluvium = alluvium,
             y = ab, label = stratum)) +
  geom_alluvium() + # aes(fill = Survived)
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal()

