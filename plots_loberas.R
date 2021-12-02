
rm(list = ls())

library(tidyverse)
library(RColorBrewer)

ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species', 'pplacer_sp') # 'Species'
sampleName <- c('Cantiles', 'Coloradito', 'Granito', 'Machos', 'Partido', 'Rasito')
dir <- '~/metagenomics/Loberas_MG/'

relationshipCol = c('P'='#6a3d9a', 'NC'='#ffff99', 'C'='#4daf4a')

featureDF <- readRDS(paste0(dir, '/featureDF.rds'))


# is_pathogenic <- function(x) (x %in% "Pathogenic")
# replaceTax <- function(x) str_replace_na(x, replacement = "Undetermined")

# featureDF %>% mutate_at(vars(ranks), funs(str_replace_na(., replacement = "Undetermined"))) -> featureDF

table(featureDF$group)

featureDF %>%
  pivot_longer(cols = ranks[-8]) %>% fill(value) %>%
  pivot_wider(names_from = name)

  # mutate_all(., funs(str_replace_na(., replacement = "Undetermined"))) -> tax

# barplot ----
facet_bar <- function(featureDF, agglom_lev = ranks[2], low_ab = 0) {
  
  featureDF %>%
    pivot_longer(all_of(sampleName), values_to = 'ab', names_to = 'Sample') %>%
    filter(ab > 0) %>%
    group_by(Sample, group) %>%
    rename( "Level" = agglom_lev) %>%
    drop_na(Level) %>%
    mutate(RA = (ab / sum(ab)) * 100) %>% # summarise(sum(RA))
    # mutate(Level = ifelse(RA <= low_ab, "ZLow", Level)) -> dataViz
    # filter(Level != 'Zlow') %>%
    filter(RA >= low_ab) %>%
    mutate(RA = (ab / sum(ab)) * 100) -> dataViz
  
  labels <- dataViz %>% pull(Level) %>% unique() %>% sort()
  colourCount = length(labels)
  
  library(ggsci)
  
  if(colourCount > 7) {
    getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
  } else
    getPalette <- pal_locuszoom()(colourCount)
  
  # if(low_ab != 0) {
  #   getPalette[length(getPalette)] <- "Black"
  #   labels[length(labels)] <- "Low abundance"
  # } 
  
  
  
  dataViz %>% group_by(group) %>%
    summarise(RA = sum(ab)) %>% 
    pivot_wider(names_from = 'group', values_from = RA) -> pout
  
  dataViz %>%
    ggplot(aes(y = RA,x = Sample, fill = Level)) +
    geom_col() +
    labs(x = "", y = "Relative Abundance (%)", 
         caption = paste0("Low Abundance (Relative Ab <", low_ab, " %)")) +
    scale_fill_manual(agglom_lev, labels = labels, values = getPalette) +
    theme_classic(base_size = 18, base_family = "GillSans") +
    facet_grid(~group) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.background = element_blank(),
          panel.border = element_blank()) -> p

  
  print(pout)
  return(p)
  
}

facet_bar(featureDF, low_ab = 1) -> p1
facet_bar(featureDF, agglom_lev = ranks[3], low_ab = 1) -> p2
facet_bar(featureDF, agglom_lev = ranks[4], low_ab = 2) -> p3

getPalette <- colorRampPalette(brewer.pal(name = 'Paired', 9))(20)

# p3 + scale_fill_manual(ranks[2], labels = labels, values = getPalette) -> p3
  

ggsave(p1, filename = "Bar_Phylum.png", path = dir, width = 10, height = 5)
ggsave(p2, filename = "Bar_Class.png", path = dir, width = 12, height = 6)
ggsave(p3, filename = "Bar_Order.png", path = dir, width = 12, height = 6)


# estimate relationship species ----
# Comensales/basales (C), asignacion no consistente (NC) y Patogena (P)

featureDF %>% group_by(group, Relationship) %>% count() %>% 
  ungroup(Relationship) %>% 
  mutate(pct = n / sum(n), pos = cumsum(pct)) %>% # summarise(sum(n))
  ggplot() + 
  aes(x = group, fill = Relationship, weight = pct) + 
  geom_bar(position = "fill") + 
  geom_label(aes(label = n, y = pos), position = position_fill(0.7), fill = "white") +
  labs(x = '', caption = 'Number of species per group in the white boxes\nInconsistent (rdp/pplacer) species assignation\nare summarised in the Undetermined group') +
  scale_fill_uchicago(name = '') -> psave

ggsave(psave, filename = 'Relationship.png', path = dir, width = 5, height = 3)

# Ocurrance upset ----

featureDF %>% filter(Relationship %in% 'Pathogenic') %>% group_by(group, Relationship) %>% count(Relationship)

featureDF %>% 
  filter(Relationship %in% 'Pathogenic') %>%
  distinct(pplacer_sp, group, Relationship) %>% 
  pivot_wider(names_from = group, values_from = group, 
              values_fn = function(x) ifelse(!is.na(x), 1, 0), values_fill = 0) -> ocurrance

ocurrance %>% count(Relationship)

cols <- names(ocurrance)[-c(1,2)]

euler <- as.data.frame(ocurrance[, -1])
rownames(euler) <- ocurrance$pplacer_sp

library(UpSetR)

n <- length(cols)
grid.col <- ggsci::pal_aaas()(n)
names(grid.col) <- cols

Myfunc <- function(row, List) {
  data <- (row["Relationship"] %in% List) 
}
# param <- list(c("Pathogenic", "Inconsistent", "non-pathogenic"))

png(filename = paste0(dir, "intersected_pathogenic_sp.png"),
    width = 780, height = 780, res = 150)

upset(euler, number.angles = 0, point.size = 3.5, line.size = 0, sets = cols, 
      keep.order = T, 
      # query.legend = "bottom",
      # nintersects = 40, 
      queries = list(
        # list(query = Myfunc, params = param, color = "orange", active = T),
        list(query = intersects, params = list(cols[1]), color = grid.col[3], active = F),
        list(query = intersects, params = list(cols[2]), color = grid.col[3], active = F),
        list(query = intersects, params = list(cols[3]), color = grid.col[3], active = F),
        list(query = intersects, params = list(cols[4]), color = grid.col[3], active = F),
        list(query = intersects, params = list(cols[5]), color = grid.col[3], active = F),
        list(query = intersects, params = list(cols[6]), color = grid.col[3], active = F)),
      mainbar.y.label = "Intersections Size", sets.x.label = "Number of Species", 
      text.scale = c(1.3, 1.5, 1, 1.5, 2, 2),
      order.by = "freq")

dev.off()

euler %>% select(cols) %>% rowSums() -> cocurrance_filtering_sp

cocurrance_filtering_sp[cocurrance_filtering_sp > 1] -> cocurrance_filtering_sp

cocurrance_filtering_sp <- names(cocurrance_filtering_sp)

writeLines(cocurrance_filtering_sp, paste0(dir, 'cocurrance_filtering_sp.list'))

featureDF %>% 
  filter(Relationship %in% 'Pathogenic') %>%
  group_by(pplacer_sp, group) %>% 
  summarise_at(vars(sampleName), sum) %>% 
  pivot_longer(cols = sampleName, names_to = 'Sample', values_to = 'ab') %>%
  group_by(pplacer_sp, group) %>% summarise(ab = sum(ab)) %>%
  pivot_wider(names_from = group, values_from = ab, values_fill = 0) %>%
  arrange(pplacer_sp) %>% filter(pplacer_sp %in% cocurrance_filtering_sp) %>% view()

# heatmap ----

featureDF %>% 
  filter(Relationship %in% 'Pathogenic') %>%
  group_by(pplacer_sp, group) %>% 
  summarise_at(vars(sampleName), sum) -> ab
# mutate(group = as.numeric(as.factor(group)))-> ab

ab %>% pivot_longer(cols = sampleName, names_to = 'Sample', values_to = 'ab') -> heatmap_df

heatmap_df %>% group_by(pplacer_sp, group) %>% summarise(ab = sum(ab)) -> heatmap_df

heatmap_df %>% pivot_wider(names_from = group, values_from = ab, values_fill = 0) -> hclust_df
hclust_df <- data.frame(hclust_df[-1], row.names = hclust_df$pplacer_sp)

hclust <- hclust(dist(hclust_df), "complete")

heatmap_df %>%
  group_by(group) %>% mutate(pct = ab/sum(ab)) %>%
  ggplot(aes(y = pplacer_sp, x = group, fill = pct)) +
  geom_raster() + 
  geom_text(aes(label = round(pct*100, digits = 2)),  vjust = 0.75, hjust = 0.5, size= 4, color = 'white') +
  # ggsci::scale_fill_gsea(name = "Membership", reverse = T, na.value = "white") +
  scale_fill_viridis_c(name = "Relative Abundance", na.value = "white", limits = c(0, 1)) +
  # ggh4x::scale_y_dendrogram(hclust = hclust, guide = ggh4x::guide_dendro(position = "left")) +
  labs(x = '', y = '') +
  guides(fill = guide_colorbar(barwidth = unit(3, "in"),
                               ticks.colour = "black", ticks.linewidth = 0.5,
                               frame.colour = "black", frame.linewidth = 0.5,
                               label.theme = element_text(size = 12))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "top") -> p1

ggsave(p1, filename = 'heatmap_pathogenic_marker.png', path = dir , width = 8, height = 10)

# heatmap ----

ab %>% pivot_longer(cols = sampleName, names_to = 'Sample', values_to = 'ab') -> heatmap_df

heatmap_df %>% pull(pplacer_sp) %>% sort(.) %>% unique() -> levels

heatmap_df %>%
  # ungroup() %>% mutate(pplacer_sp = as.factor(pplacer_sp))  %>%
  group_by(group, Sample) %>% mutate(pct = ab/sum(ab)) %>%
  arrange(pplacer_sp) %>% mutate(pplacer_sp = factor(pplacer_sp)) %>%
  # summarise(sum(pct))
  ggplot(aes(y = pplacer_sp, x = Sample, fill = pct*100)) +
  geom_raster() + 
  scale_fill_viridis_c(name = "Relative\nAbundance", na.value = "white", limits = c(0, 100)) +
  facet_grid(~ group) +
  labs(x = '', y = '') +
  guides(fill = guide_colorbar(barwidth = unit(0.3, "in"), barheight = unit(5, "in"),
                               ticks.colour = "black", ticks.linewidth = 0.5,
                               frame.colour = "black", frame.linewidth = 0.5,
                               label.theme = element_text(size = 12))) +
  theme_classic(base_size = 12, base_family = "GillSans") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank(),
        panel.border = element_blank()) -> p


ggsave(p, filename = 'heatmap_pathogenic_samples_factors.png', path = dir , width = 8, height = 10)

#

library(glue)
library(ggtext)

"<i style='color:{color}'>{pplacer_sp}</i>"

heatmap_df %>% group_by(pplacer_sp) %>% summarise(sum(ab))

heatmap_df %>%
  # distinct(pplacer_sp) %>% view() 
  group_by(group, Sample) %>%
  mutate(color = 'black') %>%
  mutate(pct = ab/sum(ab)*100, pplacer_sp = glue::glue("<i style='color:{color}'>{pplacer_sp}</i>")) %>%
  arrange(pplacer_sp) %>%
  # summarise(sum(pct))
  mutate(pct = ifelse(pct == 0, NA, pct)) %>%
  filter(pct > 0) %>%
  mutate(pplacer_sp = as.character(pplacer_sp)) %>%
  ggplot(aes(y = pplacer_sp, x = Sample, fill = pct)) +
  geom_raster() + 
  scale_y_discrete(limits=rev) +
  scale_fill_viridis_c(name = "Relative\nAbundance", 
                       na.value = "white", limits = c(0, 100)) +
  facet_grid(~ group) +
  labs(x = '', y = '') +
  guides(fill = guide_colorbar(barwidth = unit(0.3, "in"), barheight = unit(5, "in"),
                               ticks.colour = "black", ticks.linewidth = 0.5,
                               frame.colour = "black", frame.linewidth = 0.5,
                               label.theme = element_text(size = 12))) +
  theme_classic(base_size = 12) + # base_family = "GillSans"
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_markdown())  -> p


ggsave(p, filename = 'heatmap_pathogenic_samples.png', path = dir , width = 8, height = 10)

# tax tree metacoder ----

library(metacoder)
# https://grunwaldlab.github.io/metacoder_documentation/publication--05--16s_human_microbiome.html

# ab %>% pivot_wider(names_from = group, values_from = sampleName, values_fill = 0)

featureDF %>%
  filter(Relationship %in% 'Pathogenic') %>%
  select_at(vars(ranks)) %>%
  distinct(pplacer_sp, .keep_all = T) %>%
  left_join(hclust_df) %>%
  mutate_at(vars(ranks), funs(str_replace_na(., replacement = "Undetermined"))) %>%
  mutate(Species = pplacer_sp) %>% select(-pplacer_sp) %>%
  parse_tax_data(class_cols = ranks[-8], named_by_rank = TRUE) -> x

# x$data$otu_prop <- calc_obs_props(x, data = "otu_count", cols = cols)

# getting per-taxon information
x$data$tax_abund <- calc_taxon_abund(x, "tax_data", cols = cols)

x %>%
  mutate_obs("tax_abund", abundance = rowMeans(x$data$tax_abund[cols])) %>%
  filter_taxa(abundance >= 0.001, reassign_obs = c(diff_table = FALSE)) %>%
  heat_tree(node_size = n_obs,
            node_size_range = c(0.01, 0.06),
            node_size_axis_label = "Number of OTUs",
            node_color = abundance,
            node_color_axis_label = "Mean proportion of reads",
            node_label = taxon_names)

x$data$diff_table <- compare_groups(x, 
                                    data = "tax_abund", 
                                    cols = cols, 
                                    groups = cols)



set.seed(1)

output_file <- paste0(dir, "differential_heat_tree_pathogenic_sp.png")

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
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = output_file)




# test fantaxtic::fantaxtic_bar() -----
# run functions from script network_wTO_loberas.R

ps <- prepare_ps(fileNames[1], agg = T, agg_level = "Class")



library(fantaxtic)

fantaxtic_bar(ps, color_by = "Phylum", label_by = "Class", 
              palette = colorRampPalette(pal_locuszoom()(7))(10))
