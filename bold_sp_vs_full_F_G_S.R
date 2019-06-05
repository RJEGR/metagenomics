#!/usr/bin/env Rscript
# RUN FIRST bold.sp_vs_bold_full.R in order to make the alluv object:
# difine input:
# Load data ----

rm(list=ls()); 
# if(!is.null(dev.list())) dev.off()
# Clear plots

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

TL2 <- c("D", "K", "P", "C", "O", "F", "G", "S")

scale2 <- c("Domain"="#edf8b1",  "Kingdom"="#7fcdbb", "Phylum"="#2c7fb8",  
           "Class"="#feb24c",  "Order"="#addd8e",  "Family"="#31a354",
           "Genus"="#bcbddc", "Species"="#756bb1")

scale <- c("incomplete"="#E41A1C",  "complete"="#377EB8")


colnames(full.taxa.obj) <- c(TL, 'SL')
colnames(sp.taxa.obj) <- c(TL, 'SL')

# ranking labels ----

# .... 1
x_ <- NULL
for (i in 1:nrow(sp.taxa.obj)) {
  rl <- sp.taxa.obj$SL[i] + 1
  x_[[i]] <- names(sp.taxa.obj)[rl]
}

# 2.
y_ <- NULL
for (i in 1:nrow(full.taxa.obj)) {
  rl <- full.taxa.obj$SL[i] + 1
  y_[[i]] <- names(full.taxa.obj)[rl]
  
}

# 1.
x_y_rank <- data.frame(ASV = rownames(sp.taxa.obj), 
                       complete= x_, incomplete= y_, 
                       #SL_x = sp.taxa.obj$SL,
                       #SL_y = full.taxa.obj$SL,
                       x_y = sp.taxa.obj$SL - full.taxa.obj$SL,
                       stringsAsFactors = FALSE)


nrow(input_dat <- x_y_rank[x_y_rank$x_y !=0, c('ASV', 'complete','incomplete')]) # 2635

x <- melt(input_dat, id.vars  = c('ASV'),
          variable.name = 'DataBase',
          value.name = 'Rank')

levels <- c('incomplete', 'complete')
x$DataBase <- factor(x$DataBase, levels = levels)

# sanity check
colSums(aggregate(x['DataBase'], by = list(x[,'Rank']), FUN = table)[,2]) == nrow(input_dat)

# dim(out <- bbold(filter(x, Rank == TL[7:9]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 363 in filter(x, Rank == TL[2:6]
# sum(out$abund) # percent of reads with different assignation , 0.8784525 (less than 1 % of the population)
# dim(out <- out[out$abund > 1, ]) # removing singletones if no relative abundance.

# es necesario 

dim(Dom_o <- bbold(filter(x, Rank == TL[2]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 1150 asvs than reach this level in either, sp or full
dim(King_o <- bbold(filter(x, Rank == TL[3]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 1478 asvs than reach this level in either, sp or full
dim(Phyl_o <- bbold(filter(x, Rank == TL[4]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 464 asvs than reach this level in either, sp or full
dim(Class_o <- bbold(filter(x, Rank == TL[5]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 190 asvs than reach this level in either, sp or full
dim(Ord_o <- bbold(filter(x, Rank == TL[6]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 915 " "
dim(Fam_o <- bbold(filter(x, Rank == TL[7]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 109 " "
dim(Gen_o <- bbold(filter(x, Rank == TL[8]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 226 " "
dim(Sp_o <- bbold(filter(x, Rank == TL[9]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 738 " "

# nrow(Dom_o) + nrow(King_o) + nrow(Phyl_o) + nrow(Class_o) + nrow(Ord_o) + nrow(Fam_o) + nrow(Gen_o) + nrow(Sp_o)

dim(out <- rbind(data.frame(Phyl_o, Rank = TL[2]),
                 data.frame(Phyl_o, Rank = TL[3]),
                 data.frame(Phyl_o, Rank = TL[4]),
                 data.frame(Class_o, Rank = TL[5]),
                 data.frame(Ord_o, Rank = TL[6]),
                 data.frame(Fam_o, Rank = TL[7]),
                 data.frame(Gen_o, Rank = TL[8]),
                 data.frame(Sp_o, Rank = TL[9]))) # 2642 (3570 with TL[2] and TL[3])

out <- rank_n(out)

# And remove duplicates:
# ESTO TE VA A QUITAR INFORMACION EN LA CELDA RANK,
# ES NECESARIO QUEDARSE CON LA INFORMACION DEL MAX(RANK), 
# TOMAR fromLast position the duplicate element

# Ex.

tt <- out[out$ASV == 'ASV_3950',]
tt[!duplicated(tt$ASV, fromLast = TRUE) ,]

# tt %>%
#   group_by(ASV) %>%
#   filter(rank(Rank_n, ties.method = 'max')>2)

dim(out <- out[!duplicated(out$ASV, fromLast = TRUE) ,]) # 1608

table(out$Rank)

out$Rank  <- factor(out$Rank, levels = TL[2:9])

# sanity check
aggregate(out['Ref'], by = list(out[,'Rank']), FUN = table)
colSums(aggregate(out['Ref'], by = list(out[,'Rank']), FUN = table)[,2])
# complete incomplete 
# 854        754

# apply(select(out, paste("full", TL2, sep="_")), 2, n_na)
# apply(select(out, paste("sp", TL2, sep="_")), 2, n_na)

round(sum(filter(out, Ref == 'complete')$abund), digits = 2) # 1.65 (removing duplicates)
round(sum(filter(out, Ref == 'incomplete')$abund), digits = 2) # 4.45 (removing duplicates)

round(sum(out$abund)) # percent of reads within different ASV assingation between db is 6% ( with TL[2] and TL[3] )

ggplot(out, aes(seq_size, log(abund), color = Ref, shape = Ref)) +
  geom_point(alpha = 0.7, aes(size = abund)) + theme_bw() + 
  scale_color_manual(values = scale) +
  facet_wrap( ~ Rank) +
  labs(subtitle = paste0("A set of ", nrow(out), " ASVs assigned to different rank (group of x_y != 0) \n",
                         "Normalized Abundance from the total sequence abundance vs sequence size is plotted \n",
                         "The percent of ASVs belong to this set is: ",round(sum(out$abund)), "%"),
       y = 'log(RA)')

# Filtar con inner_join los resultados de out, en base a los valores de cambio x_y == 0, 
# ie. quitar del resultado out, valores que esten contenidos en la base x_y == 0 y 
# trabajar entonces con abundancia relativa

str(input_dat <- x_y_rank[x_y_rank$x_y == 0, c('ASV', 'complete','incomplete')])

# or use alluv instead of input data, alluv[, c('ASV', 'sp','full')]
y <- melt(input_dat, id.vars  = c('ASV'),
          variable.name = 'DataBase',
          value.name = 'Rank')

# 1

dim(Dom_e <- bbold(filter(y, Rank == TL[2]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 5641 asvs than reach this level in either, sp or full
dim(King_e <- bbold(filter(y, Rank == TL[3]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 4445 asvs than reach this level in either, sp or full
dim(Phyl_e <- bbold(filter(y, Rank == TL[4]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 523 asvs than reach this level in either, sp or full
dim(Class_e <- bbold(filter(y, Rank == TL[5]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 196 asvs than reach this level in either, sp or full
dim(Ord_e <- bbold(filter(y, Rank == TL[6]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 449 " "
dim(Fam_e <- bbold(filter(y, Rank == TL[7]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 169 " "
dim(Gen_e <- bbold(filter(y, Rank == TL[8]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 316 " "
dim(Sp_e <- bbold(filter(y, Rank == TL[9]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 6904 " "
# 2
dim(equals <- rbind(data.frame(Dom_e, Rank = TL[2]),
                    data.frame(King_e, Rank = TL[3]),
                    data.frame(Phyl_e, Rank = TL[4]),
                    data.frame(Class_e, Rank = TL[5]),
                    data.frame(Ord_e, Rank = TL[6]),
                    data.frame(Fam_e, Rank = TL[7]),
                    data.frame(Gen_e, Rank = TL[8]),
                    data.frame(Sp_e, Rank = TL[9]))) # 8557 (18643 with TL[2] and TL[3])

dim(equals <- equals[!duplicated(equals$ASV, fromLast = TRUE) ,]) # just in case

round(sum(equals$abund)) # percent of reads within equal ASV assingation between db is 84 %

equals$Rank <- factor(equals$Rank, levels = TL[2:9])
# # sanity check ()
round(sum(equals$abund)) + round(sum(out$abund)) + 10 == 100


# visualize seq_size and abundance :
equals_sbt <- select(equals, ASV, Ref, Rank, abund, seq_size)
equals_sbt$Ref <- "Both"

nrow(seq_size_vs_abund <- rbind(equals_sbt, select(out, ASV, Ref, Rank, abund, seq_size)))

# sanity check below fails
# nrow(full.taxa.obj) == nrow(sp.taxa.obj) #TRUE
# nrow(out) + nrow(equals)
# (nrow(equals) + nrow(out)) - nrow(sp.taxa.obj) # 1027
# check missed ASVs
nrow(missed <- subset(x_y_rank, !(ASV %in% seq_size_vs_abund$ASV)))

missed <- select(missed, -x_y)

dim(missed_tbl <- bbold(missed, fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 1027 " "

round(sum(missed_tbl$abund))

missed_tbl$Rank <- missed_tbl$Ref

missed_sbt <- data.frame(select(missed_tbl, ASV, Ref, Rank, abund, seq_size))

missed_sbt$Ref <- "Missed"

dim(missed_sbt <- missed_sbt[!duplicated(missed_sbt$ASV, fromLast = TRUE) ,])

nrow(seq_size_vs_abund <- rbind(seq_size_vs_abund, missed_sbt))
nrow(sp.taxa.obj)

# compare 
ggplot(seq_size_vs_abund, aes(seq_size, log(abund), color = Rank)) +
  geom_point(alpha = 0.5, aes(size = abund)) + theme_bw(base_size = 12) + 
  scale_color_manual(scale = scale) +
  facet_wrap( ~ Ref) +
  scale_color_manual(values=c(scale2, scale)) + coord_flip() +
  labs(caption = paste0("Here's present a set of ", nrow(out), " and ", nrow(equals), " different ASVs assigned to R rank level due to sequence resolution [x_y != 0] and [x_y == 0], respectively\n",
                         "Normalized Abundance from the total sequence abundance vs sequence size is plotted \n",
                         "The percent of ASVs belong to this set is: ",round(sum(out$abund)), "% and ", round(sum(equals$abund)), " %, respectively",
                          "A set of ", nrow(missed), " ASVs were missed from the analysis"))

# 1.
test0 <- melt(out, id.vars = c('ASV', 'seq_size', 'abund', 'x_y', 'Ref', 'Rank'), variable.name = 'DataBase', value.name = 'lineage')
test0 <- select(test0, ASV, Ref, Rank, lineage)

# 2.
test1 <- melt(equals, id.vars = c('ASV', 'seq_size', 'abund', 'x_y','Ref', 'Rank'), variable.name = 'DataBase', value.name = 'lineage')
test1 <- select(test1, ASV, Ref, Rank, lineage)

nrow(test1 <- test1[!duplicated(test1$lineage),]) # 1130 unique lineage 

# test1$DataBase <- factor(test1$DataBase, levels = paste("sp", TL2[4:9], sep="_"))

# Revisamos que las asignaciones difertes sean realmente nuevas entre bases y no se encuentren
# contenidas en el objeto equals

test0 %>% 
  group_by(lineage) %>%
  anti_join(test1, by = "lineage")  %>%
  #select(ASV, lineage, Rank) %>%
  as.data.frame() -> test

# dim(test) # 374 (422 with TL[2] and TL[3])

dim(test <- test[!duplicated(test$ASV),]) # 142

# sanity check
# should print data.frame of 0 columns, or TRUE in == 0
# due to non anti_join lineage from out (test0) in equals (test1) object
ncol(test1[test$ASV %in% test1$ASV]) == 0

# In comparison from the full-set of 2635 [table(alluv$Ref)], 1083 (complete) and 1552 (incomplete)
# only a set of 142 ASVs (47 [complete] and 95 [incomplete]) shows different and unique assigment due to database composition 
table(test$Ref) < table(alluv$Ref)
colSums(aggregate(test[,'Ref'], by=list(test[,'Rank']), FUN = table)[2])

library(tidyverse)
# head(spread(test0, Rank, Ref))
# head(spread(test1, Rank, Ref))

upset.obj <- test
#dim(upset.obj0 <- spread(upset.obj, Ref, Rank))

upset.obj[upset.obj$Ref == 'incomplete', 'Ref'] <- 2
upset.obj[upset.obj$Ref == 'complete', 'Ref'] <- 1
dim(upset.obj0 <- spread(upset.obj, Rank, Ref))
upset.obj0[is.na(upset.obj0)] <- 0
upset.obj0 <- upset.obj0[,-c(1,2)]

library(UpSetR)
# venn diagrams
# https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html

# movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
#                    header = T, sep = ";")

#upset(upset.obj0, number.angles = 30)
# 
# upset(movies, nsets = 7, number.angles = 30, point.size = 3.5, line.size = 2, 
#       mainbar.y.label = "Genre Intersections", sets.x.label = "Movies Per Genre", 
#       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))

dim(out0 <- bbold(test, fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 142

# resumen: ----
# de 2635 asvs con asignacion diferente `nrow(x_y_rank[x_y_rank$x_y !=0, c('ASV', 'sp','full')])`
# solo 1130 `nrow(test1 <- test1[!duplicated(test1$lineage),])`, fueron ASVs con diferencias entre las bases
# la diferencia de 2635-1130 corresponde a ASVs con linaje contenido en alguna de las dos bases
# apesar de demostrar asignacion diferente, por tanto no representan informacion perdida. 
# Mientras los 1130 representan informacion ganda en alguna de las dos bases 
# 


# if(ncol(test1[test$ASV %in% test1$ASV]) == 0) {
#   out %>%
#   group_by(ASV) %>%
#   select( -Rank) %>%
#   inner_join(test, by = 'ASV') %>%
#   select( -Rank) %>%
#   as.data.frame() -> out0
# } else {out0 <- out}


# 0. Select Family
data0 <- rbind(
  data.frame(table(out0$complete_F), db = 'complete'),
  data.frame(table(out0$incomplete_F), db = 'incomplete')
  
)

names(data0) <- c('lineage', 'Size', 'DataBase')
data0$Rank <- TL[7]

# 1. Select Genus
data1 <- rbind(
  data.frame(table(out$complete_G), db = 'complete'),
  data.frame(table(out$incomplete_G), db = 'incomplete')
  
)

names(data1) <- c('lineage', 'Size', 'DataBase')
data1$Rank <- TL[8]

# 2, Select Order

data2 <- rbind(
  data.frame(table(out$complete_S), db = 'complete'),
  data.frame(table(out$incomplete_S), db = 'incomplete')
  
)

names(data2) <- c('lineage', 'Size', 'DataBase')
data2$Rank <- TL[9]

# 3. Parse results
dim(ntaxa_data <- rbind(data0,data1, data2)) # 254
ntaxa_data$Rank <- factor(ntaxa_data$Rank, levels = TL[7:9])
ntaxa_data <- ntaxa_data[order(-ntaxa_data$Size),]

# aggregate colsums by 


tax_sp <- select(as.data.frame(out0), abund, paste("complete", TL2, sep="_"))
tax_full <- select(as.data.frame(out0), abund, paste("incomplete", TL2, sep="_"))

# Calculate abundance ----
# coherence with data from ntaxa size
# by

Family <- rbind( data.frame(aglom_ab(tax_sp, 'complete_F'),
                            DataBase = 'complete',
                            Rank = 'Family'),
                 data.frame(aglom_ab(tax_full, 'incomplete_F'),
                            DataBase = 'incomplete',
                            Rank = 'Family'))
# By 
Genus <- rbind( data.frame(aglom_ab(tax_sp, 'complete_G'),
                           DataBase = 'complete',
                           Rank = 'Genus'),
                data.frame(aglom_ab(tax_full, 'incomplete_G'),
                           DataBase = 'incomplete',
                           Rank = 'Genus'))
# By 

Species <- rbind( data.frame(aglom_ab(tax_sp, 'complete_S'),
                           DataBase = 'complete',
                           Rank = 'Species'),
                data.frame(aglom_ab(tax_full, 'incomplete_S'),
                           DataBase = 'incomplete',
                           Rank = 'Species'))
# Plot aglomerated abundance ----

abund_data <- rbind(Family, Genus, Species)
abund_data$Rank <- factor(abund_data$Rank, levels = TL[7:9])
abund_data <- abund_data[order(-abund_data$Size),]


# merge ntaxa and abund
ntaxa_data$abund <- 'nASVs'
abund_data$abund <- "nreads_pct"

data <- rbind(ntaxa_data, abund_data)

data$abund <- factor(data$abund, levels = c('nreads_pct', 'nASVs'))
data$DataBase <- factor(data$DataBase, levels = levels)

ggplot(filter(data, abund == 'nASVs'),  aes(lineage, Size, fill = DataBase, group = DataBase)) +
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") + coord_flip() +
  geom_text(aes(label = Size, color = DataBase), size = 3, position = position_dodge(width = 0.5)) +
  scale_color_brewer(palette = "Set1") +
  #facet_grid(Rank~ abund, scales = 'free') +
  facet_wrap(~Rank, scales = 'free') +
  theme_classic(base_size=7) +
  labs(y = 'n ASVs',
    title = 'Number of taxa during database assignation')
    #caption = paste0('Using the subset: filter(x, Rank == "Domain","Kingdom","Phylum","Class", "Order") \n Singletones removed'))

# 2. only species 

ggplot(filter(data, Rank == 'Species'),  aes(lineage, Size, fill = DataBase, group = DataBase)) +
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") + coord_flip() +
  # geom_text(aes(label = Size, color = DataBase), size = 3, position = position_dodge(width = 0.5)) +
  geom_text(data = subset(data, Rank == 'Species' && abund == 'nASVs'), aes(label=Size,  color = DataBase), 
                position = position_dodge(width = 0.5), check_overlap =FALSE, size = 3) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(Rank~ abund, scales = 'free') +
  theme_bw(base_size=10) +
  labs(y = 'n Size',
       title = 'Number of taxa during database assignation',
       subtitle = paste0('Using the subset: filter(data, Rank == "Species", abund == "nreads_pct") \n', 
                        'The percent of Relative abundance is ', 
                        round(sum(filter(data, Rank == 'Species', abund == 'nreads_pct')$Size), digits = 2), '%'))


# search some rare species assignations
lineage <- filter(data, abund == 'nASVs' & DataBase == 'complete' & Rank == 'Species')$lineage

length(lineage)

out0[out0$incomplete_S %in% lineage, paste0('complete_', TL2) ]
# df1[with(df1, grepl("B|F", paste(Col2, Col3))),]

# calcular radio de abundancia de ntaxa:nread ----
# 161:17438/161

count.tbl0 <- read.table(paste0(path_BOLD, "/",count_tbl))
# se tienen que contar incluso singletones, sum(count.tbl0[rowSums(count.tbl0) > 1,])

# ternary plot
# see figure 3 from https://www.biorxiv.org/content/biorxiv/early/2019/03/12/575928.full.pdf to replicate it
# install.packages('ggtern')

# install.packages('Ternary')
#library('Ternary')
# TernaryPlot() ----

library(ggtern)

tern_obj <- subset(x_y_rank_m, x_y != 0)
tern_obj0 <- aggregate(tern_obj[,'DataBase'], by=list(tern_obj[,'Rank']), FUN = table)

x_y_rank$complete <- factor(x_y_rank$complete, levels = TL[-1])
tern_shared <- data.frame(table(c(select(subset(x_y_rank, x_y == 0), complete))))

datavis <- data.frame(Rank = tern_obj0[,1], tern_obj0[,2], shared = tern_shared$Freq)
datavis <- data.frame(Rank = tern_obj0[,1], apply(datavis[-1], 2, function(x) { x / sum(x) * 100})) 

my_theme <- function(base_size = 12, base_family = "") {
  theme_custom(base_size, base_family, col.T = '#377EB8', col.L = '#E41A1C', 
               col.R = "darkgreen") + theme(tern.panel.background = element_rect(fill = "white"), 
                                            tern.panel.grid.minor = element_line(color = "gray90"), 
                                            tern.axis.arrow.show = TRUE)
}

ggtern(data=datavis,aes(incomplete, complete, shared)) + 
  #geom_point() + 
  my_theme() +
  geom_text(aes(label = Rank)) +
  scale_color_manual(scale = scale)


# okay, tenemos dos posibles hipotesis:
# - a medida que la base de datos aumenta, el algorimo desconoce que hacer con tantas secuencias 
# - a medida que la base de datos aumenta, el algorimo enriquece la informacion para hacer una asignacion mas confiable
# - a medida que la base de datos disminuye, el algoritmo introduce falsos positivos pues asigna la secuencia mas parecida,
# - a medida que la base de datos disminuye, el algoritmo introduce falsos negativos pues no logra asignar el taxon 'unclassify'

# trabajar con el set de secuencias de los asvs unicos (~ 50) y comparar su asignacion entre las dos bases de referencia:
# - evaluar distancia hamming
# - un blast entre las bases de datos y la lista de asvs
# - 
# remover asvs en base a tamano de secuencia,
# codon de paro (a.a), 
# distancia

# seleccionamos asvs ----
# track the list of asvs from databases (full, sp), and use sequences to analyse hamming distance or other
# this task (or other) is neccesary to report false positive or false negative (also define this analysis thorug sanger mock68 )

asv_subset <- out$ASV

seqs0 <- readDNAStringSet(paste0(path_BOLD,'/',fasta_file))
seqs <- seqs0[names(seqs0) %in% asv_subset]
seqs <- seqs[match(asv_subset, names(seqs)),]

# to rename, use last lineage assigned

x_ <- llineage(sp.taxa.obj)
y_ <- llineage(full.taxa.obj)

x_y_lineage <- data.frame(ASV = rownames(sp.taxa.obj), sp= x_, full= y_)
x_y_l_subset <- x_y_lineage[x_y_lineage$ASV %in% asv_subset, ]
x_y_l_subset <- x_y_l_subset[match(asv_subset, x_y_l_subset$ASV),]


library(ShortRead)

# seqs_headers <- paste0(">filtFs_Read", seq(length(seqs)))

#file_name <- "blastn_nr/blastnr_mergers.fasta"
#save_fasta <- c(rbind(seqs_headers, dfM$sequence))
#write(save_fasta, file=file_name)

# dfM$id <- seqs_headers

# for hamming
write(out$ASV, file=paste0(path_BOLD,"/non_shared_asvs.out" ))



taxonomy.file <- paste0(path_BOLD,'/',bold_all)
taxonomy.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
taxonomy <- as.data.frame(t(taxonomy))
tax0 <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)

rownames(tax0) <- taxonomy.obj$V1

tax <- unite(tax0, col = 'lineage', sep = ';')

tax_subset <- data.frame(ASV = rownames(tax)[which(rownames(tax) %in% asv_subset)], 
                         lineage = tax[which(rownames(tax) %in% asv_subset), ])

write(tax_subset$lineage, file=paste0(path_BOLD,"/non_shared_taxa.out" ))

#  get the taxonomy, then de db id and sequence from db
# awk 'NR==FNR{a["ASV"$0];next}/^ASV/{f=0;}($0 in a)||f{print;f=1}' non_shared_asvs.out run014_t2_ASVs.BOLD_public_species.wang.taxonomy > test
#  awk '{print $1}' non_shared_asvs.out | xargs -I {} grep "{}" run014_t2_ASVs.BOLD_public_species.wang.taxonomy
# awk 'NR==FNR{a[$1]; next} {for (i in a) if (index($0, i)) print $1}' non_shared_asvs.out run014_t2_ASVs.BOLD_public_species.wang.taxonomy
#### BLAST by hand against nt. Defaults. Megablast.
#### Save hit table (text) (using F2PFFZXG015-Alignment-HitTable.csv)
# and load functions from mock_analysis_hamming.R

seqsM <- as.data.frame(seqs)


allM <- matrix(0, ncol=nrow(seqsM), nrow=1)

colnames(allM) <- rownames(seqsM)
rownames(allM) <- 'abund'

allM['abund',] <- out$abund

dfM <- data.frame(t(allM))

dfM$sequence <- seqsM$x

# comprare vs ref ----
# search in taxonomy.*.csv
dir = '/Users/cigom/metagenomics/db/bold/BOLD_public_trim'
fasta.file = 'dada2_asv/run012_20190329_COI/mock_hits_relax_ASVs.fasta'
reference.file = 'ictio_coi_sanger114.fasta'

ref <- readFasta(reference.file)
refSeqs <- as.character(sread(ref))
strain <- as.character(id(ref))

dfM.vs.ref <- outer(dfM$sequence, refSeqs, evalDist, band=-1)
dfM$refdist <- apply(dfM.vs.ref, 1, min)

# continue with

blast_file <- '/Users/cigom/metagenomics/COI/run012/mock_P68/mock_parameter_definition/blastn_nr/blastnr_mergers-Alignment.txt'
blast_file <- '/Users/cigom/metagenomics/COI/species_resolution_per_db/F2PFFZXG015-Alignment.txt'

dfM$hit <- isHit100(dfM, blast_file)
dfM$oo <- isOneOff(dfM, blast_file)


