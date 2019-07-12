# 
# Julio 2019
# Ricardo Gomez

rm(list = ls())

#install.packages("bold")
require('bold')
citation(package = 'bold')

dim(res <- bold_specimens(taxon='Cumacea')) # heavy data

head(res[,1:8])

dim(res <- bold_specimens(taxon='Cumacea'))

# large data
library("taxize")

x <- downstream("Arthropoda", db = "ncbi", downto = "class")

nms <- x$Arthropoda$childtaxa_name

checks <- bold_tax_name(nms)

# all is good
checks[,1:5]
out <- lapply(nms, bold_seq)

# or from
# PUBLIC DATA PORTAL - BIN LIST
# http://boldsystems.org/index.php/Public_BINSearch?query=Mexico

dir <- '/Users/cigom/metagenomics/COI/BOLD_SYSTEM_MEX/'

tbl <- read_delim(paste0(dir, 'bold_data_mex.txt'), delim = '\t')

# filter only COI records
# c('COII', 'COXIII', 'COI-3P')
COI_5P <- filter(tbl, markercode == 'COI-5P')
COI_3P <- filter(tbl, markercode == 'COI-3P')
COXIII <- filter(tbl, markercode == 'COXIII')
COII <- filter(tbl, markercode == 'COII')

nrow(tbl_COI <- rbind(COI_5P, COI_3P, COXIII, COII)) # 29285

# get field of interest:

taxonomy <- c('bin_uri', 'phylum_name', 'class_name', 'order_name', 'family_name', 'genus_name', 'species_name')

location <- c('lat', 'lon', 'country', 'region')

seq_properties <- c('sequenceID', 'markercode', 'genbank_accession', 'nucleotides', 'directions', 'seq_primers')

# 
# Taxonomy ----
# 

tbl_select <- select(tbl_COI, taxonomy)

# bins 
bins <- data.frame(table(tbl_select$bin_uri))

head(bins <- bins[order(bins$Freq, decreasing = TRUE),])

# bin_test <- filter(tbl_select, bin_uri == 'BOLD:ACO9738')
# 
# table(bin_test$phylum_name)
phylum_name <- data.frame(table(tbl_select$phylum_name))
class_name <- data.frame(table(tbl_select$class_name))
order_name <- data.frame(table(tbl_select$order_name))

# head(class_name <- class_name[order(class_name$Freq, decreasing = TRUE),])

#
# location ----
#

tbl_select <- select(tbl_COI, location, seq_properties) # taxonomy

coors <- read.csv("~/Downloads/metadata_xiximi06.csv", header=TRUE, stringsAsFactors = FALSE)
coors$Transect <- substr(coors$Station, 1,1)
cicese <- data.frame(lat = coors$Latitude, lon = coors$Longitude,
                     country = 'Mexico', region = 'GoM')
# lat    lon country region


table(tbl_select$country)
# 29091
library(ggrepel)
library(ggplot2)
library(tidyverse)

mex <- ggplot2::map_data("world", region = "Mexico") 

usa <- map_data("state") %>%
  subset(region %in% c("florida", "texas", "louisiana", "mississippi", "alabama"))

world <- rbind(mex, usa)

gg <- ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill = NA, color = "#D2D2D2", size = 0.7) +  # 
  # coord_fixed(xlim = c(-100,-80), ylim = c(15,30)) + # GoM only
  coord_fixed(xlim = c(-120,-80), ylim = c(15,31)) + # GoM only
  theme_bw()

# filter(world, long < -95 & long < -90 & lat < 20 & lat > 15)

gg2 <- gg + 
  # geom_point(data = tbl_select, aes(lon, lat), alpha = 0.7, color = "#365B60") +
  geom_point(data = tbl_select, aes(lon, lat, color = markercode), alpha = 0.7) +
  scale_color_brewer(palette = "RdYlBu") +
  #geom_text_repel(data = coors,
  #               aes(Longitude, Latitude, label = Station),
  #                segment.color = "grey50") + 
  theme(text = element_text(size=12)) +
  theme_bw() +
  ggtitle("")

gg2 + geom_point(data = cicese, aes(lon, lat), alpha = 0.7, shape = 2, color = '#365B60')

# 
# Sequence ----
tbl_select$seq_len <- nchar(tbl_select$nucleotides)

ggplot(tbl_select, aes(x=seq_len)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(seq_len)),
             color="blue", linetype="dashed", size=0.5) +
  theme_classic()

p <- ggplot(tbl_select, aes(x=seq_len, color = markercode)) + 
  geom_density() + scale_color_brewer(palette="Dark2") +
  facet_wrap(~ markercode, scales = 'free') +
  theme_classic()

## search database ncbi for:
tbl_select <- select(tbl_COI, seq_properties)

# with genbank accession number
table(is.na(tbl_select$genbank_accession))
# FALSE  TRUE 
# 11,456 17,829

markercode <- data.frame(table(tbl_select$markercode))
head(markercode <- markercode[order(markercode$Freq, decreasing = TRUE),])

# Var1  Freq
# 2 COI-5P 29091
# 1 COI-3P   190
# 3   COII     2
# 4 COXIII     2

# 


