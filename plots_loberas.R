
rm(list = ls())

library(tidyverse)

ranks <- c('Kingdom',  'Phylum',  'Class',  'Order', 'Family', 'Genus', 'Species', 'pplacer_sp') # 'Species'
sampleName <- c('Cantiles', 'Coloradito', 'Granito', 'Machos', 'Partido', 'Rasito')
dir <- '~/metagenomics/Loberas_MG/'

relationshipCol = c('P'='#6a3d9a', 'NC'='#ffff99', 'C'='#4daf4a')

featureDF <- readRDS(paste0(dir, '/featureDF.rds'))

featureDF %>% 
  mutate(Relationship = recode_factor(Relationship, P = 'Pathogenic', 
                                      NC = 'Inconsistent', C = 'Commensal')) -> featureDF

# is_pathogenic <- function(x) (x %in% "Pathogenic")
# replaceTax <- function(x) str_replace_na(x, replacement = "Undetermined")

# featureDF %>% mutate_at(vars(ranks), funs(str_replace_na(., replacement = "Undetermined"))) -> featureDF



featureDF %>%
  pivot_longer(cols = ranks[-8]) %>% fill(value) %>%
  pivot_wider(names_from = name) %>% view()

table(featureDF$Relationship)
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
    mutate(Level = ifelse(RA <= low_ab, "ZLow", Level)) -> dataViz
  
  labels <- dataViz %>% pull(Level) %>% unique() %>% sort()
  colourCount = length(labels)
  
  library(ggsci)
  
  if(colourCount > 7) {
    getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
  } else
    getPalette <- pal_locuszoom()(colourCount)
  
  if(low_ab != 0) {
    getPalette[length(getPalette)] <- "Black"
    labels[length(labels)] <- "Low abundance"
  } 
  
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
          panel.border = element_blank())
  
}

facet_bar(featureDF, low_ab = 0) -> p

ggsave(p, filename = "Bar_Phylum.png", path = dir, width = 10, height = 5)

# estimate relationship species ----
# Comensales/basales (C), asignacion no consistente (NC) y Patogena (P)

featureDF %>% group_by(group, Relationship) %>% count() %>% 
  ungroup(Relationship) %>% 
  mutate(pct = n / sum(n), pos = cumsum(pct)) %>% summarise(sum(n))
  ggplot() + 
  aes(x = group, fill = Relationship, weight = pct) + 
  geom_bar(position = "fill") + 
  geom_label(aes(label = n, y = pos), position = position_fill(.7), fill = "white") +
  labs(x = '', caption = 'Number of species per group in the white boxes\nInconsistent (rdp/pplacer) species assignation\nare summarised in the Undetermined group') +
  scale_fill_uchicago(name = '') -> psave

ggsave(psave, filename = 'Relationship.png', path = dir, width = 5, height = 3)

# Ocurrance upset ----

featureDF %>% 
  filter(Relationship %in% 'Pathogenic') %>%
  distinct(pplacer_sp, group, Relationship) %>% # count(Relationship)
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
param <- list(c("Pathogenic", "Undetermined", "Commensal"))

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
        list(query = intersects, params = list(cols[4]), color = grid.col[3], active = F)),
      mainbar.y.label = "Intersections Size", sets.x.label = "Number of Species", 
      text.scale = c(1.3, 1.5, 1, 2, 2, 2),
      order.by = "freq")

dev.off()



# tax tree metacoder ----

library(metacoder)

featureDF %>% 
  filter(Relationship %in% 'Pathogenic') %>%
  group_by(pplacer_sp, group) %>% 
  summarise_at(vars(sampleName), sum) -> ab
  # mutate(group = as.numeric(as.factor(group)))-> ab

# ab %>% pivot_wider(names_from = group, values_from = sampleName, values_fill = 0)

featureDF %>%
  filter(Relationship %in% 'Pathogenic') %>%
  select_at(vars(ranks)) %>%
  distinct(pplacer_sp, .keep_all = T) %>%
  left_join(ocurrance) %>%
  mutate_at(vars(ranks), funs(str_replace_na(., replacement = "Undetermined"))) %>%
  mutate(Species = pplacer_sp) %>% select(-pplacer_sp) %>%
  parse_tax_data(class_cols = ranks[-8], named_by_rank = TRUE) -> x

# getting per-taxon information
x$data$tax_abund <- calc_taxon_abund(x, "tax_data", cols = cols)

x$data$diff_table <- compare_groups(x, 
                                    data = "tax_abund", 
                                    cols = cols, 
                                    groups = cols)

set.seed(1)

output_file <- paste0(dir, "differential_heat_tree.pdf")

heat_tree_matrix(x,
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



