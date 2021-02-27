rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

files <- list.files(path = dir, pattern = ".xlsx$", full.names = TRUE)
mtd_file <- list.files(path = dir, pattern = "CIAD_MappingFile.txt$", full.names = TRUE)

library(readxl)
library(tidyverse)
library(ggsci)

obj <- read_xlsx(files[2])

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
  mutate_at(ranks, 
            funs(str_replace_all(., c("_[1-9]" = "")))) -> obj

obj %>%
  select_if(is.double) %>%
  names() -> colNames

barTax <- function(obj, colNames, agglom_lev = "Phylum", low_ab = 1) {
  obj %>% 
    pivot_longer(cols = colNames, names_to = "Index", values_to = 'ab') %>%
    filter(ab > 0) %>%
    inner_join(mtd) %>%
    rename( "Level" = agglom_lev) %>%
    group_by(Level, Tissue) %>%
    summarise(ab = sum(ab), Freq = length(Level > 0)) %>%
    group_by(Tissue) %>%
    mutate(RA = (ab / sum(ab)) * 100) %>%
    mutate(Level = ifelse(RA < low_ab, "ZLow", Level)) -> dataViz
  
  labels <- dataViz %>% pull(Level) %>% unique() %>% sort()
  colourCount = length(labels)
  
  if(colourCount > 7) {
    getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
  } else
    getPalette <- pal_locuszoom()(colourCount)
  
  
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

PPhylum <- barTax(obj, colNames, "Phylum")
PClass <- barTax(obj, colNames, "Class", low_ab = 1)
POrder <- barTax(obj, colNames, "Order", low_ab = 1)

ggsave(PPhylum, filename = "Bar_Phylum.png", path = dir, 
       width = 10, height = 8)
ggsave(PClass, filename = "Bar_Class.png", path = dir, 
       width = 10, height = 8)
ggsave(POrder, filename = "Bar_Order.png", path = dir, 
       width = 10, height = 8)

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
