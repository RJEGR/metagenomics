#!/usr/bin/env Rscript

# rm(list=ls()); 
# Clear plots

if(!is.null(dev.list())) dev.off()

# Functions ----

options(stringsAsFactors = FALSE)
#

source(file = "~/Documents/GitHub/metagenomics/readtx.R")

library(ggplot2)
library(reshape2)
library(dplyr)


# Set filenames ----
path_BOLD <- '/Users/cigom/metagenomics/COI/species_resolution_per_db'
bold_all <- 'run014_t2_ASVs.ALL.wang.taxonomy'
bold_sp <- 'run014_t2_ASVs.BOLD_public_species.wang.taxonomy'
fasta_file <- 'run014_t2_ASVs.fasta'
count_tbl <- 'run014_t2_ASVs_count.table'
# Load files ----

full.taxa.obj <- read_rdp(paste0(path_BOLD,'/',bold_all))
sp.taxa.obj <- read_rdp(paste0(path_BOLD,'/',bold_sp))


TL <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

colnames(full.taxa.obj) <- c(TL, 'SL')
colnames(sp.taxa.obj) <- c(TL, 'SL')



# SL Level Resolution vis ----

SL <- rbind(data.frame(SL = full.taxa.obj$SL,
                 DataBase = 'incomplete'),
            data.frame(SL = sp.taxa.obj$SL,
                       DataBase = 'complete'))


# SL <- SL[SL$SL > 0 ,]
SL$SL <- factor(SL$SL, levels = names(table(SL$SL)))

# and
SL_agg <- aggregate(SL[,'DataBase'], by=list(SL[,'SL']), FUN = table)
SL_agg <- data.frame(level = SL_agg[,1], SL_agg[,2])


# sanity check
asv_size <- nrow(sp.taxa.obj)
if(nrow(full.taxa.obj) == asv_size) {
  if(colSums(SL_agg[-1])[1] == asv_size)
    if(colSums(SL_agg[-1])[2] == asv_size)
    SL_agg$Rank <- TL[-1]
    
  SL_aggM <- melt(SL_agg, id.vars = )
  SL_aggM$Rank <- factor(SL_aggM$Rank, levels = TL[-1])
}

SL_aggM$pct <- round(c((SL_agg$incomplete / asv_size ) * 100, (SL_agg$complete / asv_size ) * 100), digits = 0)


levels <- c("incomplete","complete")

SL_aggM$variable <- factor(SL_aggM$variable, levels = levels)


# taxonomuc resolution ----

ggplot(SL_aggM, aes(x=Rank, y=value, group=variable, fill=variable, color=variable)) + 
  geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") +
  geom_text(data = SL_aggM, aes(label=pct), size=7, vjust=1, color = 'black', check_overlap = TRUE) +
  scale_color_brewer(palette = "Set1") +
  labs(title="ASVs fully and insufficiently identified due to database composition",
       subtitle = paste0('The size of the ASV sample is: ', nrow(sp.taxa.obj) ),
       x = 'Taxonomy level',
       y = '# ASVs') +
  theme_bw(base_size = 14)
       
# non_sp selection ----


# if(identical(rownames(sp.taxa.obj), rownames(full.taxa.obj))) {
#   non_sp <- rownames(sp.taxa.obj[sp.taxa.obj$SL != max(sp.taxa.obj$SL),])
#   sp_non_sp <- as_tibble(sp.taxa.obj[rownames(sp.taxa.obj) %in% non_sp, ])
#   full_non_sp <- as_tibble(full.taxa.obj[rownames(full.taxa.obj) %in% non_sp, ])
# }

#  quantitative comparison ----

n_na <- function(x) ( sum(is.na(x))) # sum(is.na(sp_non_sp[,2]))

na_number <- data.frame(
  complete = apply(sp.taxa.obj[ -10], 2, n_na),
  incomplete = apply(full.taxa.obj[-10], 2, n_na))

na_number$Level <- rownames(na_number)
na_number$diff <- na_number$incomplete - na_number$complete
na_number_m <- melt(na_number, id.vars = c('Level', 'diff'),
               variable.name = 'DataBase',
               value.name = 'n')

na_number_m$Level <- factor(na_number_m$Level, levels = TL)
na_number_m$DataBase <- factor(na_number_m$DataBase, levels = levels)

ggplot(subset(na_number_m, diff != 0), aes(Level, n, fill = DataBase, group = DataBase)) + 
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  geom_text(data = subset(na_number_m, DataBase == 'complete' & diff != 0), aes(label=abs(diff)), size=4, vjust=1, color = 'black') +
  scale_fill_brewer(palette = "Set1") +
  labs(x="Rank", y="# ASVs (Undetermined | Unclassified)", 
       title = 'Number of NAs ASVs across fully and infufficiently identified BOLD databases') +
  theme_bw(base_size = 12)

# count the level of change ----
x <- sp.taxa.obj$SL
y <- full.taxa.obj$SL

diff_x_y <- data.frame(
            Change = -6:6,# data.frame(table(x - y))[,1],
            Freq = data.frame(table(x - y))[,2],
            stringsAsFactors = FALSE)
# 

diff_x_y$Ref <- 'Shared'
diff_x_y[diff_x_y[,1] < 0, 'Ref'] <- 'incomplete'
diff_x_y[diff_x_y[,1] > 0, 'Ref'] <- 'complete'


diff_x_y$Change <- abs(diff_x_y$Change)
diff_x_y$Change <- factor(diff_x_y$Change, levels = 1:6)
diff_x_y$Ref <- factor(diff_x_y$Ref, levels = levels)

diff_val <- sum(diff_x_y$Freq) -  max(diff_x_y$Freq)

ggplot(subset(diff_x_y, Change != 0), aes(x=Change, y=log10(Freq), group=Ref, fill=Ref, color=Ref)) + 
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") +
  geom_text(data = subset(diff_x_y, Change != 0), aes(label=Freq), size=4, vjust=1, color = 'black') +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(title= paste0("Number of ASVs with difference in the assignation"),
       subtitle = paste0("The number of shared ASVs are: ", max(diff_x_y$Freq), " and ASVs with change in the assignation due to DataBase: ",diff_val) ) +
  theme_bw(base_size = 12)

# or
#if(full_non_sp[,1] == 'root') {full_non_sp <- full_non_sp[-1]}
#if(sp_non_sp[,1] == 'root') {sp_non_sp <- sp_non_sp[-1]}

# ranking labels ----

# .... 1
x_ <- lrank(sp.taxa.obj)

# 2.
y_ <- lrank(full.taxa.obj)

# 1.
x_y_rank <- data.frame(ASV = rownames(sp.taxa.obj), 
                            complete= x_, incomplete= y_, 
                       #SL_x = sp.taxa.obj$SL,
                       #SL_y = full.taxa.obj$SL,
                       x_y = sp.taxa.obj$SL - full.taxa.obj$SL,
                       stringsAsFactors = FALSE)

# sanity check of shared asvs
table(select(filter(x_y_rank, x_y == 0), complete)) == table(select(filter(x_y_rank, x_y == 0), incomplete))

# 2.
x_y_rank_m <-  melt(filter(x_y_rank, x_y != 0), id.vars = c('ASV', 'x_y'),
                     variable.name = 'DataBase',
                     value.name = 'Rank')

x_y_rank_m$Rank <-  factor(x_y_rank_m$Rank, levels = TL[-1])


x_y_rank_m$DataBase <-  factor(x_y_rank_m$DataBase, levels = levels)

# Level change by rank
# agrupa los niveles por Rank y entonces ploteas el panel_grid

ggplot(x_y_rank_m, 
       aes(y=..count.., x=abs(x_y), fill = DataBase)) + 
  geom_histogram(aes(y=..count..), position=position_dodge(), alpha=0.5, bins = 10) +
  #geom_density(alpha=.7) + 
  scale_x_continuous(name = "Change", breaks  =  c(1:6)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x="change", y="Frequency of ASVs", 
       title = 'Level of Resolution between full and sp  BOLD dabases') +
  facet_wrap( ~ Rank , scales = 'free_y') +
theme_bw(base_size = 12)

# and the zero across Ranks
    
# alluvial

dim(alluv <- subset(x_y_rank, x_y != 0)) # 2635 ASVs with difference in assignation over databases
alluv <- db_color(alluv) 
table(alluv$Ref) # 1083 (complete) and 1552 (incomplete)

# alluv <- alluv[alluv$ASV %in% sample(alluv$ASV, 20) ,]

alluv$complete <-  factor(alluv$complete, levels = TL[-1])
alluv$incomplete <-  factor(alluv$incomplete, levels = TL[-1]) 

alluv <- with(alluv, alluv[order(complete),])
alluv <- select(alluv, complete, incomplete, Ref)

# RColorBrewer::brewer.pal(3, "Set1")

scale <- c("incomplete"="#E41A1C",  "complete"="#377EB8")

library(ggalluvial)

alluv_long <- to_lodes_form(data.frame(alluv), key = "x", axes = 1:2)

# Sanity check
colSums(aggregate(alluv_long['x'], by = list(alluv_long[,'stratum']), FUN = table)[,2]) == nrow(alluv)


is_alluvia_form(alluv_long, silent = TRUE)

ggplot(data = alluv_long,
       aes(x = x, stratum = stratum, alluvium = alluvium,
           label = stratum)) +
  geom_alluvium(aes(fill = Ref)) +
  geom_stratum() + geom_text(stat = "stratum") +
  theme_minimal() + 
  scale_fill_manual(values=scale) +
  ggtitle("The frequency distribution over Ranks.") +
  xlab("Database") + ylab("Number of ASVs")

# https://github.com/corybrunson/ggalluvial/issues/13
# For some strata appear at multiple aces warning messages

# (Previous version)
# ggplot(data = alluv,
#        aes(axis1 = complete, axis2 = incomplete, color = Ref)) +
#   stat_stratum(geom = "errorbar", na.rm = TRUE, alpha = 0.5) +
#   geom_line(stat = "alluvium") +
#   stat_alluvium(geom = "pointrange") +
#   scale_color_manual(values=scale) +
#   geom_text(stat = "stratum", label.strata = TRUE, color = 'black') +
#   scale_x_discrete(limits = c("complete", "incomplete")) +
#   theme_minimal()

# 1) 


scale <- c("Domain"="#edf8b1",  "Kingdom"="#7fcdbb", "Phylum"="#2c7fb8",  
           "Class"="#feb24c",  "Order"="#addd8e",  "Family"="#31a354",
           "Genus"="#bcbddc", "Species"="#756bb1")

levels <- c('complete', 'incomplete')
alluv_long$Ref <- factor(alluv_long$Ref, levels = levels) 

ggplot(data = alluv_long,
       aes(x = x, stratum = stratum, alluvium = alluvium,
           label = stratum)) +
  geom_stratum() + geom_text(stat = "stratum") +
  geom_flow(aes(fill = stratum, alpha = Ref), stat = "alluvium",
            aes.bind = FALSE, lode.guidance = "rightward") +
  theme_minimal() + scale_fill_manual(values=scale) +
  ggtitle("The frequency distribution over Ranks.") +
  xlab("Database") + ylab("Number of ASVs")

# (Previous version)
# ggplot(data = alluv,
#        aes(axis1 = complete, axis2 = incomplete)) + 
#   geom_alluvium() + 
#   geom_stratum(na.rm = TRUE) + 
#   scale_x_discrete() +
#   geom_text(stat = "stratum", label.strata = TRUE, size = 3) +
#   scale_x_discrete(limits = c("complete", "incomplete"), expand = c(.05, .05)) +
#   # theme(legend.position = "bottom") +
#   theme_minimal() + 
#   # use of lode controls
#   geom_flow(aes(fill = complete, alpha = Ref), stat = "alluvium",
#             aes.bind = FALSE, lode.guidance = "rightward") +
#             # color = "black",  linetype="dashed") +
#   scale_fill_manual(values=scale)

# and 2

levels <- c('incomplete', 'complete')

alluv_long$Ref <- factor(alluv_long$Ref, levels = levels) 

ggplot(data = alluv_long,
       aes(x = x, stratum = stratum, alluvium = alluvium,
           label = stratum)) +
  geom_stratum() + geom_text(stat = "stratum") +
  geom_flow(aes(fill = stratum, alpha = Ref), stat = "alluvium",
            aes.bind = FALSE, lode.guidance = "rightward") +
  theme_minimal() + scale_fill_manual(values=scale) +
  ggtitle("The frequency distribution over Ranks.") +
  xlab("Database") + ylab("Number of ASVs")

# filtering rare assingation ----
# Analsis 2 ----
# siguiente script (bold_sp_full_F_G_S.R)

# A) ventajas de la full (incomplete)
# cuando se procesa la base de datos species se implementa el filtro en r complete.cases, 
# ie. la celda debe contener valores distintos a NA; bold, contiene individuos con asignacion
# hasta especie con 'gaps' entre los niveles taxonomicos, por lo que este tipo de asignaciones son
# desconocidas en la base de datos de bold-species Ej: reconocemos que hay asignaciones hasta especie
# nuevas en la base de datos bold-full y que en bold-sp se asignaron por debajo de este nivel:

# filter(x, full == 'Species' & x_y == -6)
# dim(out <- bbold(filter(x, full == 'Species' & x_y == -6), fasta_file = fasta_file, count_tbl = count_tbl))
# out[out$full_C == 'Polychaeta',] # <- example

# B) ventaja de la spp
# Un mayor numero de ASVs son asignados a especie usando la sp, mientras que en full, 
# se quedan entre los niveles Orden en adelante,
# existe la hipotesis de que la complejidad de la base de datos full ocasiona una serie de falsos positivos
# es necesario conocer cual es el limite al cual el algoritmo tiene su maximo rendimiento (rdp)


# EXISTEN NUEVOS INDIVIDUOS USANDO LA BASE DE DATOS FULL?
# length(full_subset <- x[x$full == TL[2:6], 'ASV'])






# exit ----

# Relative abundances:
count.tbl0 <- read.table(paste0(path_BOLD, "/",count_tbl))
count.obj <- data.frame(apply(count.tbl0, 2, function(x) (x / sum(x)) * 100 ))
# sanity check

table(colSums(count.obj) == 100) # 49 samples

y <- x$ASV
subset <- y[!duplicated(y)]

count_vis <- count.obj[rownames(count.obj) %in% full_subset, ]

tax <- select(out, paste("sp", TL2, sep="_"))
# tax <- select(out, paste("full", TL2, sep="_"))

library(phyloseq)

phyloseq = phyloseq(otu_table(count_vis, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')))

physeq <- prune_taxa(taxa_sums(phyloseq) > 0.1, phyloseq) # # Removing any taxa of abundance of 1

# keepTaxa = apply(X = as(otu_table(phyloseq), "matrix") > 0, # Removing any abundance of zero
#                 MARGIN = 1, FUN = sum) > 1 # Remove ASVs not k greater than k (2L)
# clean.obj = prune_taxa(keepTaxa, phyloseq)

Family <-tax_glom(physeq, taxrank="sp_F")
Family = subset_taxa(Family, sp_F!="NA")

library(superheat)


heat.tbl <- otu_table(Family)
rownames(heat.tbl) <- tax_table(Family)[,'sp_F']
# colnames(heat.tbl) <- sample_data(Family)$`Estación`

rowSums(heat.tbl)


# apply(heat.tbl, 1, function(x) (x  * sum(x) ))
library(superheat)
superheat(heat.tbl, 
          # scale the variables/columns
          scale = FALSE,
          # change the color
          heat.pal = c( "#FFFFD9", "#081D58"),
          # change the color of the labels
          left.label.col = "white",
          # change the size of the label text
          left.label.text.size = 3,
          # change the angle of the label text
          bottom.label.text.angle = 90,
          bottom.label.text.size = 3,
          bottom.label.text.alignment = "right",
          bottom.label.col = "white",
          # add barplot next to the rows
          # Titles
          column.title = 'Muestras',
          column.title.size = 4,
          row.title = 'Taxones',
          row.title.size = 4,
          # remove the grid
          grid.hline = FALSE,
          grid.vline = FALSE,
          # retain original order of rows/cols
          pretty.order.rows = TRUE,
          #pretty.order.cols = TRUE,
          col.dendrogram  = TRUE
)

# ggplot heatmap
ggplot(data = xplot,
       aes(x = variable, y = Genes, fill = value)) + 
  geom_raster(hjust = 0, vjust = 0) +
  labs(fill = "log2(FPKM)", x = 'Samples', y ='Genes') +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1), axis.text.y = element_text(size = 5))

# #
# plot_bar(Family, x= "Sample", fill = 'sp_F') +
#   geom_bar(stat = "identity", position = "stack") +
#   # scale_fill_manual(values = getPalette(colourCount)) +
#   ylab("Abundancia") + coord_flip() +
#   guides(fill=guide_legend(ncol=3)) +
#   theme(panel.background = element_blank()) +
#   labs(x = 'Estación', caption = "")

