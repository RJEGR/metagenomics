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
source(file = "~/Documents/GitHub/metagenomics/hamming.R")

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
x_ <- lrank(sp.taxa.obj)
y_ <- lrank(full.taxa.obj)

# 2.
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

nrow(Dom_o) + nrow(King_o) + nrow(Phyl_o) + nrow(Class_o) + nrow(Ord_o) + nrow(Fam_o) + nrow(Gen_o) + nrow(Sp_o)

dim(out <- rbind(data.frame(Dom_o, Rank = TL[2]),
                 data.frame(King_o, Rank = TL[3]),
                 data.frame(Phyl_o, Rank = TL[4]),
                 data.frame(Class_o, Rank = TL[5]),
                 data.frame(Ord_o, Rank = TL[6]),
                 data.frame(Fam_o, Rank = TL[7]),
                 data.frame(Gen_o, Rank = TL[8]),
                 data.frame(Sp_o, Rank = TL[9]))) # 5270 (3570 with TL[2] and TL[3])

out$Rank <- factor(out$Rank, levels = TL[-1])
out$Ref <- factor(out$Ref, levels = levels)
# out <- rank_n(out)

# And remove duplicates:
# ESTO TE VA A QUITAR INFORMACION EN LA CELDA RANK,
# ES NECESARIO QUEDARSE CON LA INFORMACION DEL MAX(RANK), 
# TOMAR fromLast position the duplicate element

# Ex.
head(out$ASV[which(table(out$ASV) >= 2)])
tt <- out[out$ASV == 'ASV_17358', c('ASV', 'Ref','Rank')]
tt[!duplicated(tt$ASV, fromLast = TRUE) ,]

# keep the last ASV with highes assignation

dim(out <- out[!duplicated(out$ASV, fromLast = TRUE) ,]) # 2635

table(out$Rank)


# sanity check
aggregate(out['Ref'], by = list(out[,'Rank']), FUN = table)

table(out$Ref)

# complete incomplete 
# 1552        1083

# apply(select(out, paste("full", TL2, sep="_")), 2, n_na)
# apply(select(out, paste("sp", TL2, sep="_")), 2, n_na)

round(sum(filter(out, Ref == 'complete')$abund), digits = 2) # 3.08 (removing duplicates)
round(sum(filter(out, Ref == 'incomplete')$abund), digits = 2) # 11.04 (removing duplicates)

round(sum(out$abund)) # percent of reads within different ASV assingation between db is 14% ( with TL[2] and TL[3] )

ggplot(out, aes(seq_size, -log(abund), color = Ref, shape = Ref)) +
  geom_point(alpha = 0.7, aes(size = abund)) + theme_bw() + 
  scale_color_manual(values = scale) +
  facet_wrap( ~ Rank) +
  labs(subtitle = paste0("A set of ", nrow(out), " ASVs assigned to different rank (group of x_y != 0) \n",
                         "Normalized Abundance from the total sequence abundance vs sequence size is plotted \n",
                         "The percent of ASVs belong to this set is: ",round(sum(out$abund)), "% of the RA"),
       y = '-log(RA %)')

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
sum(equals$abund) + sum(out$abund)  == 100


# visualize seq_size and abundance :
equals_sbt <- select(equals, ASV, Ref, Rank, abund, seq_size)
equals_sbt$Ref <- "Both"

nrow(seq_size_vs_abund <- rbind(equals_sbt, select(out, ASV, Ref, Rank, abund, seq_size))) # 20251 ASVs

# sanity check
nrow(full.taxa.obj) == nrow(sp.taxa.obj) #TRUE
nrow(out) + nrow(equals) == nrow(sp.taxa.obj)

# if missed, check missed ASVs
#nrow(missed <- subset(x_y_rank, !(ASV %in% seq_size_vs_abund$ASV))) # 1027
#table(missed$incomplete)
# missed <- select(missed, -x_y)
# round(sum(missed_tbl$abund)) # 7
# missed_tbl$Rank <- missed_tbl$Ref
# missed_sbt <- data.frame(select(missed_tbl, ASV, Ref, Rank, abund, seq_size))
# missed_sbt$Ref <- "Missed"
# dim(missed_sbt <- missed_sbt[!duplicated(missed_sbt$ASV, fromLast = TRUE) ,])

# compare 
ggplot(seq_size_vs_abund, aes(seq_size, -log(abund), color = Rank)) +
  geom_point(alpha = 0.5, aes(size = abund)) + theme_bw(base_size = 12) + 
  scale_color_manual(scale = scale) +
  facet_wrap( ~ Ref) +
  scale_color_manual(values=scale2) + coord_flip() +
  labs(caption = paste0("Here's present a set of ", nrow(out), " and ", nrow(equals), " different ASVs assigned to R rank level due to sequence resolution [x_y != 0] and [x_y == 0], respectively\n",
                         "Normalized Abundance from the total sequence abundance vs sequence size is plotted \n",
                         "The percent of ASVs belong to this set is: ",round(sum(out$abund)), "% and ", round(sum(equals$abund)), " % of the RA, respectively"))

# 1. out - change
test0 <- melt(out, id.vars = c('ASV', 'seq_size', 'abund', 'x_y', 'Ref', 'Rank'), variable.name = 'DataBase', value.name = 'lineage')
test0 <- select(test0, ASV, Ref, Rank, lineage)

# 2. equals
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


head(test$ASV[which(table(test$ASV) >= 3)])
tt <- x[x$ASV == 'ASV_8000',]
tt[!duplicated(tt$ASV, fromLast = TRUE) ,]

dim(test <- test[!duplicated(test$ASV, fromLast = TRUE),]) # 142

aggregate(test['Ref'], by = list(test[,'Rank']), FUN = table)

# sanity check
# should print data.frame of 0 columns, or TRUE in == 0
# due to non anti_join lineage from out (test0) in equals (test1) object

ncol(test1[test$ASV %in% test1$ASV]) == 0

# In comparison from the full-set of 2635 [table(alluv$Ref)], 1083 (complete) and 1552 (incomplete)
# only a set of 142 ASVs (47 [complete] and 95 [incomplete]) shows different and unique assigment due to database composition 
sum(table(test$Ref)) # < sum(table(alluv$Ref))

colSums(aggregate(test[,'Ref'], by=list(test[,'Rank']), FUN = table)[2])

# library(tidyverse)
# head(spread(test0, Rank, Ref))
# head(spread(test1, Rank, Ref))
# 
# upset.obj <- test
# #dim(upset.obj0 <- spread(upset.obj, Ref, Rank))

# upset.obj[upset.obj$Ref == 'incomplete', 'Ref'] <- 2
# upset.obj[upset.obj$Ref == 'complete', 'Ref'] <- 1
# dim(upset.obj0 <- spread(upset.obj, Rank, Ref))
# upset.obj0[is.na(upset.obj0)] <- 0
# upset.obj0 <- upset.obj0[,-c(1,2)]

# library(UpSetR)
# venn diagrams
# https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html

# movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
#                    header = T, sep = ";")

#upset(upset.obj0, number.angles = 30)
# 
# upset(movies, nsets = 7, number.angles = 30, point.size = 3.5, line.size = 2, 
#       mainbar.y.label = "Genre Intersections", sets.x.label = "Movies Per Genre", 
#       text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))

dim(out0_c <- bbold(filter(test, Ref == 'complete'), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 142
dim(out0_i <- bbold(filter(test, Ref == 'incomplete'), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 142

#round(sum(out0$abund), digits = 2)  # 1.26 % of reads
round(sum(out0_c$abund), digits = 2)  # 0.37 % of reads
round(sum(out0_i$abund), digits = 2)  # 0.89 % of reads


# if(ncol(test1[test$ASV %in% test1$ASV]) == 0) {
#   out %>%
#   group_by(ASV) %>%
#   select( -Rank) %>%
#   inner_join(test, by = 'ASV') %>%
#   select( -Rank) %>%
#   as.data.frame() -> out0
# } else {out0 <- out}

# Produce ASV size object
# Select Order


# dim(out0_c <- filter(out0, Ref == 'complete'))
# dim(out0_i <- filter(out0, Ref == 'incomplete'))

data0 <- rbind(
  data.frame(table(out0_c$complete_O), db = 'complete'),
  data.frame(table(out0_i$incomplete_O), db = 'incomplete')
  
)

# 0. Select Family
data1 <- rbind(
  data.frame(table(out0_c$complete_F), db = 'complete'),
  data.frame(table(out0_i$incomplete_F), db = 'incomplete')
  
)
# 1. Select Genus
data2 <- rbind(
  data.frame(table(out0_c$complete_G), db = 'complete'),
  data.frame(table(out0_i$incomplete_G), db = 'incomplete')
  
)


# 2, Select Species

data3 <- rbind(
  data.frame(table(out0_c$complete_S), db = 'complete'),
  data.frame(table(out0_i$incomplete_S), db = 'incomplete')
  
)

names(data0) <- c('lineage', 'Size', 'DataBase')
names(data1) <- c('lineage', 'Size', 'DataBase')
names(data2) <- c('lineage', 'Size', 'DataBase')
names(data3) <- c('lineage', 'Size', 'DataBase')

data0$Rank <- TL[6]
data1$Rank <- TL[7]
data2$Rank <- TL[8]
data3$Rank <- TL[9]



# 3. Parse results
dim(lineage_tbl <- rbind(data0,data1, data2, data3)) # 126 different lineages
lineage_tbl$Rank <- factor(lineage_tbl$Rank, levels = TL[6:9])
#lineage_tbl$DataBase <- factor(lineage_tbl$Ref, levels = levels)
lineage_tbl <- lineage_tbl[order(-lineage_tbl$Size),]


# Produce Relative abundance object
# aggregate colsums by 

dim(tax_sp <- select(as.data.frame(out0_c), abund, paste("complete", TL2, sep="_")))
dim(tax_full <- select(as.data.frame(out0_i), abund, paste("incomplete", TL2, sep="_")))

# Calculate abundance ----
# coherence with data from ntaxa size


Order <- rbind( data.frame(aglom_ab(tax_sp, 'complete_O'),
                            DataBase = 'complete',
                            Rank = 'Order'),
                 data.frame(aglom_ab(tax_full, 'incomplete_O'),
                            DataBase = 'incomplete',
                            Rank = 'Order'))
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

abund_tbl <- rbind(Order, Family, Genus, Species)
abund_tbl$Rank <- factor(abund_tbl$Rank, levels = TL[6:9])
abund_tbl <- abund_tbl[order(-abund_tbl$Size),]


aggregate(test['Ref'], by = list(test[,'Rank']), FUN = table)
aggregate(lineage_tbl['DataBase'], by = list(lineage_tbl[,'Rank']), FUN = table)
colSums(aggregate(abund_tbl['DataBase'], by = list(abund_tbl[,'Rank']), FUN = table)[2])


# merge ntaxa and abund
lineage_tbl$abund <- 'nASVs'
abund_tbl$abund <- "nreads_pct"

dim(data <- rbind(lineage_tbl, abund_tbl)) # 252

data$abund <- factor(data$abund, levels = c('nreads_pct', 'nASVs'))
data$DataBase <- factor(data$DataBase, levels = levels)

ggplot(data,  aes(lineage, Size, fill = DataBase, group = DataBase)) +
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") + coord_flip() +
  #geom_text(aes(label = Size, color = DataBase), size = 3, position = position_dodge(width = 0.5)) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(Rank~ abund, scales = 'free') +
  # facet_wrap(~Rank, scales = 'free') +
  theme_bw(base_size=7) +
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

# calcular radio de abundancia de ntaxa:nread ----
# 161:17438/161

# count.tbl0 <- read.table(paste0(path_BOLD, "/",count_tbl))
# se tienen que contar incluso singletones, sum(count.tbl0[rowSums(count.tbl0) > 1,])

# ternary plot
# see figure 3 from https://www.biorxiv.org/content/biorxiv/early/2019/03/12/575928.full.pdf to replicate it
# install.packages('ggtern')

# install.packages('Ternary')
#library('Ternary')
# TernaryPlot() ----

library(ggtern)


tern_obj <-  melt(filter(x_y_rank, x_y != 0), id.vars = c('ASV', 'x_y'),
                    variable.name = 'DataBase',
                    value.name = 'Rank')

table(out$Ref)
nrow(out_c <- filter(out, Ref == 'complete'))
nrow(out_i <- filter(out, Ref == 'incomplete'))

out_c$Rank <- factor(out_c$Rank, levels = TL[-1])
out_i$Rank <- factor(out_i$Rank, levels = TL[-1])
# 
# tern_obj <- subset(x_y_rank_m, x_y != 0)
# tern_obj0 <- aggregate(tern_obj[,'DataBase'], by=list(tern_obj[,'Rank']), FUN = table)

tern_c <- data.frame(table(out_c$Rank))
tern_i <- data.frame(table(out_i$Rank))

x_y_rank$complete <- factor(x_y_rank$complete, levels = TL[-1])
tern_shared <- data.frame(table(c(select(subset(x_y_rank, x_y == 0), complete))))

datavis <- data.frame(Rank = tern_c$Var1, complete = tern_c$Freq, incomplete = tern_i$Freq, shared = tern_shared$Freq)

#datavis <- data.frame(Rank = tern_obj0[,1], tern_obj0[,2], shared = tern_shared$Freq)
datavis <- data.frame(Rank = tern_c[,1], apply(datavis[-1], 2, function(x) { x / sum(x) * 100})) 

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

# search some rare species assignations
lineage <- filter(data, abund == 'nASVs' & DataBase == 'complete' & Rank == 'Species')$lineage

length(lineage) # 20

# out0[out0$complete_S %in% lineage, paste0('complete_', TL2) ]
# df1[with(df1, grepl("B|F", paste(Col2, Col3))),]


# track the list of asvs from databases (full, sp), and use sequences to analyse hamming distance or other
# this task (or other) is neccesary to report false positive or false negative (also define this analysis thorug sanger mock68 )

dim(out0 <- bbold(test, fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = FALSE)) # 142
asv_subset <- out0$ASV
length(asv_subset <- asv_subset[!duplicated(asv_subset)]) # 142

seqs0 <- readDNAStringSet(paste0(path_BOLD,'/',fasta_file))
seqs <- seqs0[names(seqs0) %in% asv_subset]
seqs <- seqs[match(asv_subset, names(seqs)),]

# to rename, use last lineage assigned

x_ <- llineage(sp.taxa.obj)
y_ <- llineage(full.taxa.obj)

dim(x_y_lineage <- data.frame(ASV = rownames(sp.taxa.obj), complete= x_, incomplete= y_))
dim(x_y_l_subset <- x_y_lineage[x_y_lineage$ASV %in% asv_subset, ])
x_y_l_subset <- x_y_l_subset[match(asv_subset, x_y_l_subset$ASV),]
# seq(length(seqs))
if(identical(names(seqs), x_y_l_subset$ASV))
  if(identical(names(seqs), out0$ASV)){
  seqs_headers <- paste0(">", names(seqs),"::", x_y_l_subset$complete, "::", x_y_l_subset$incomplete, "::", out0$Ref)
  seqs_headers <- str_replace_all(seqs_headers, c(" "="_"))
  
  save_fasta <- c(rbind(seqs_headers, as.data.frame(seqs)$x))
}


file_name <- paste0(path_BOLD, "/", "non_shared_asvs.fasta")
write(save_fasta, file=file_name)


# for hamming
# write(out$ASV, file=paste0(path_BOLD,"/non_shared_asvs.out" ))

taxonomy.file <- paste0(path_BOLD,'/',bold_all)
taxonomy.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
taxonomy <- as.data.frame(t(taxonomy))
tax0 <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)

tax0 <- mutate_all(data.frame(tax0), funs(str_replace_all(., c("_unclassified"="", "Unclassified"=""))))

rownames(tax0) <- taxonomy.obj$V1

tax <- tidyr::unite(tax0, col = 'lineage', sep = ';')

tax_subset <- data.frame(ASV = rownames(tax)[which(rownames(tax) %in% asv_subset)], 
                         lineage = tax[which(rownames(tax) %in% asv_subset), ])

nrow(tax_subset <- tax_subset[match(asv_subset, tax_subset$ASV),])
non_shared_out <- tax_subset$lineage

non_shared_out <- non_shared_out[!duplicated(non_shared_out)]

write(non_shared_out, file=paste0(path_BOLD,"/non_shared_asvs.taxonomy" ))

#  get the taxonomy, then de db id and sequence from db
# awk 'NR==FNR{a["ASV"$0];next}/^ASV/{f=0;}($0 in a)||f{print;f=1}' non_shared_asvs.out run014_t2_ASVs.BOLD_public_species.wang.taxonomy > test
#  awk '{print $1}' non_shared_asvs.out | xargs -I {} grep "{}" run014_t2_ASVs.BOLD_public_species.wang.taxonomy
# awk 'NR==FNR{a[$1]; next} {for (i in a) if (index($0, i)) print $1}' non_shared_asvs.out run014_t2_ASVs.BOLD_public_species.wang.taxonomy
#### BLAST by hand against nt. Defaults. Megablast.
#### Save hit table (text) (using F2PFFZXG015-Alignment-HitTable.csv)
# and load functions from mock_analysis_hamming.R

#dfM is a data.frame of ASV-sequence with abundance from both, complete and incomplete list

ids <- paste0(names(seqs),"::", x_y_l_subset$complete, "::", x_y_l_subset$incomplete, "::", out0$Ref)

ids <- str_replace_all(ids, c(" "="_"))

nrow(dfM <- data.frame(id = ids, ASV = names(seqs), complete = x_y_l_subset$complete, incomplete = x_y_l_subset$incomplete, 
                  Ref = out0$Ref, abund = out0$abund, seq_size = out0$abund, sequence = as.data.frame(seqs)$x))


# str(lineage_out <- c(x_y_l_subset$complete, x_y_l_subset$incomplete))
# select unique names of lineage than differ between a databse vs database

str(lineage_out <- c(x_y_l_subset$complete, x_y_l_subset$incomplete))
str(lineage_out <- lineage_out[!duplicated(lineage_out)]) # 72 lineages

# search in db*.csv
dir = '/Users/cigom/metagenomics/db/bold/BOLD_public_trim/'
fasta.file = 'non_shared_asvs.fasta'
reference.file = 'non_shared_asvs.vs.BOLD_public.ALL.fasta'

library(ShortRead)

# make the subset of BOLD database
# grep -F -f non_shared_asvs.taxonomy BOLD_public.ALL.tax > non_shared_asvs.vs.db.taxonomy &
# awk '{print $2}' non_shared_asvs.vs.db.taxonomy | sort | uniq > non_shared_asvs.vs.db.located.taxonomy
# diff -u non_shared_asvs.vs.db.located.taxonomy non_shared_asvs.taxonomy | grep -E "^\+" | wc -l
# awk '{print $1}' non_shared_asvs.vs.db.taxonomy > non_shared_asvs.vs.db.taxonomy.ids
# awk 'NR==FNR{a[">"$0];next}/^>/{f=0;}($0 in a)||f{print;f=1}' non_shared_asvs.vs.db.taxonomy.ids BOLD_public.ALL.fasta > non_shared_asvs.vs.BOLD_public.ALL.fasta

ref <- readFasta(paste0(dir, reference.file))

ref <- readDNAStringSet(paste0(dir, reference.file))

strain <- as.character(names(ref))

refSeqs0 <- data.frame(ref)$ref

# head(alphabetFrequency(ref, baseOnly = FALSE))
# refSeqs <- str_replace_all(refSeqs0, c("-"="", "[+]"="", "[.]"="",
#                                       "M"="", "R"="", "W"="",
#                                       "S"="", "Y"="", "K"="", "V"="", "H"="",
#                                       "D"="", "B"="", "N"="")) # remove gaps

# A   C   G   T M R W S Y K V H D B N   - + .
# [1,] 152 170 105 247 0 0 0 0 0 0 0 0 0 0 1 226 0 0
# [2,] 139 101 107 185 0 0 0 0 0 0 0 0 0 0 0 126 0 0

alphab <- data.frame(id=strain, alphabetFrequency(ref, baseOnly=TRUE))

length(refSeqs_id <- alphab[alphab$other ==0, 'id'])

refSeqs0 <- ref[names(ref) %in% refSeqs_id]
refSeqs0 <- refSeqs0[match(refSeqs_id, names(refSeqs0)),]

head(alphabetFrequency(refSeqs0, baseOnly=TRUE))

refSeqs <- data.frame(refSeqs0)$refSeqs0

save(list = ls(all=TRUE), file = paste0(path_BOLD, "hamming.RData"))
# load(paste0(path_BOLD, "hamming.RData"))
# load function evalDist
# source(file = "~/Documents/GitHub/metagenomics/hamming.R")
dfM.vs.ref <- outer(dfM$sequence, refSeqs, evalDist, band=-1)

dfM$refdist <- apply(dfM.vs.ref, 1, min)

save(list = ls(all=TRUE), file = paste0(path_BOLD, "hamming.RData"))


hist(log10(dfM$refdist))

sum(dfM$refdist==0 & dfM$Ref == 'complete') # 3
sum(dfM$refdist==0 & dfM$Ref == 'incomplete') # 2

# blast parsing # <---- we're here!!!

blast_file <- paste0(path_BOLD, '/', 'FJRA2ND7014-Alignment.txt')


dfM$hit <- isHit100(dfM, blast_file)
dfM$oo <- isOneOff(dfM, blast_file)

# 1e-10
get_ham_nln(data.frame(sequence=dfM$sequence[dfM$E_10>0], abundance=dfM$E_10[dfM$E_10>0]))

dfM$abund <- out0$abund
dfM$seq_size = out0$seq_size

dfM$ham_nln <- as.integer(NA)

dfM$ham_nln[dfM$abund>0] <- get_ham_nln(data.frame(sequence=dfM$sequence, abundance=dfM$abund[dfM$abund>0]))

# define accuracy
dfM$Accuracy <- getAccuracy(dfM)


# Visualize  -----

scl.y <- scale_y_log10(limits=c(10^-6, 10^-0.5))
scl.x <- scale_x_log10(limits=c(1, 90))
thm <- theme(plot.margin=rep(unit(0, "in"),4), panel.grid=element_blank(), legend.key=element_blank())
accScale <- scale_shape_manual(name="Accuracy", values=c("Reference"=0, "Exact"=2, "One Off"=4, "Other"=8))
lab.x <- xlab("Hamming (Log-scale)")
lab.y <- ylab("Frequency (Log-scale)")

#
totM.ex <- sum(dfM$abund)

m1 <- ggplot(data=dfM[dfM$abund>0,], aes(x=ham_nln, y=ham_nln/totM.ex, shape=Accuracy))
m1 <- m1 + geom_point(size=3)
m1 <- m1 + accScale + theme_bw() + ggtitle("")
m1 <- m1 + thm + lab.x + lab.y
# m1 <- m1 + scl.y + scl.x + thm + lab.x + lab.y
m1

acc_tbl <- table(dfM$Accuracy[dfM$abund>0])


xaScale <- scale_color_manual(name="Vs. ", values=c("Same"="black", "Added"="#0099FF", "Lost"="grey70"))
dfM$xa <- "Same"
dfM$xa[dfM$abund==0 & dfM$abund>0] <- "Added"
dfM$xa[dfM$E_120>0 & dfM$E_10==0] <- "Lost"
dfM$xa[dfM$E_120==0 & dfM$E_10==0] <- "N/A"

m2 <- ggplot(data=dfM[dfM$abund>0 | dfM$abund>0,], aes(x=ham_nln, y=E_10/totM.ex, shape=Accuracy, color=xa))
m2 <- m2 + geom_point(size=3)
m2 <- m2 + geom_point(data=dfM[dfM$E_120>0 & dfM$E_10==0,], aes(x=E_120, y=E_120/totM.ex), size=2, show.legend=FALSE)
m2 <- m2 + accScale + xaScale + theme_bw() + ggtitle("OMEGA 1E-120")
m2 <- m2 + scl.y + scl.x + thm + lab.x + lab.y
m2 <- m2 + geom_vline(xintercept=0.03*nchar(dfM$sequence[[1]]), linetype="dashed", color="grey25")
m2
