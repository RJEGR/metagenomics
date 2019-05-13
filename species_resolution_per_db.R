#!/usr/bin/env Rscript

rm(list=ls()); 

# Load libraries ----
.cran_packages <- c("stringr", "reshape2", "tidyverse") 
# .cran_packages <- c("stringr", "tidyverse", "reshape2", "ggplot2") 

# 1. cran_package
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}

# Load packages into session, and print package version
sapply(c(.cran_packages), require, character.only = TRUE)


args = commandArgs(trailingOnly=TRUE)

# wangtax = args[1]
# tag <- strsplit(wangtax, "[.]")[[1]][1]
# wang tax from:
path = '/Users/cigom/metagenomics/COI/species_resolution_per_db'

bold.file = "run014_t2_ASVs.BOLD_public_species_v02.wang.taxonomy"
midori.file  = "run014_t2_ASVs.midori_unique_DB_0.wang.taxonomy"
coarbitrator.file = "run014_t2_ASVs.Coarbitrator_COI_nuc_curated.wang.taxonomy"
# To load file:  wang.taxonomy
source(file = "~/Documents/GitHub/metagenomics/converted2worms/cigom_load_dada2_tables.R")
setwd(path)
## Load assignation result ----

bold <- read.wangtax(bold.file)
midori <- read.wangtax(midori.file)
coarbitrator <- read.wangtax(coarbitrator.file)

TL <- names(midori)[-c(1,ncol(midori))]

if (identical(names(bold), names(midori) )) { TL <- names(bold)[-c(1,ncol(bold))]
    } 

# Contamos el numero de indeterminados por base de datos a nivel filo ----
n_undetermined <- function(x) ( nrow(x[x[,3] == 'Undetermined',]) )

print(n_undetermined_bold <- n_undetermined(bold))
print(n_undetermined_midori <- n_undetermined(midori))
print(n_undetermined_coarbitrator <- n_undetermined(coarbitrator))


# Back the SL position per database ----
coarbitrator_ <- NULL

for (i in 1:nrow(coarbitrator)) {
  rl <- coarbitrator$SL[i] + 1
  coarbitrator_[[i]] <- list(rank=names(coarbitrator)[rl], linage=coarbitrator[i,rl]) }

coarbitrator_ <- do.call(rbind, coarbitrator_)

midori_ <- NULL

for (i in 1:nrow(midori)) {
    rl <- midori$SL[i] + 1
    midori_[[i]] <- list(rank=names(midori)[rl], linage=midori[i,rl]) }

midori_ <- do.call(rbind, midori_)
# 
bold_ <-NULL
for (i in 1:nrow(bold)) {
    rl <- bold$SL[i] + 1
    bold_[[i]] <- list(rank=names(bold)[rl], linage=bold[i,rl]) }

bold_ <- do.call(rbind, bold_)


# compare the last parent between db ----

names(midori)[-1] <- paste("midori", names(midori)[-1], sep="_")
names(bold)[-1] <- paste("bold", names(bold)[-1], sep="_")
names(coarbitrator)[-1] <- paste("coarbitrator", names(coarbitrator)[-1], sep="_")

# bold vs midori ----

LCR <- data.frame(midori, bold[,-1], coarbitrator[,-1], diff_x_y =  midori[,9] - bold[,9], diff_x_z = coarbitrator[,9] - bold[,9], stringsAsFactors =  FALSE)
# aniadir una columna, rango 1
LCR$bold_vs_midori <- 1
# recorrer hacia el 7 en loop
for(t in 1:7){
    # Generar nombre de la columna con taxonomia por morfologia
    midori_rank <- paste("midori", TL[t], sep="_")
    # Generar nombre de la columna con taxonomica por molecular
    bold_rank <- paste("bold", TL[t], sep="_")
    # Aquellos renglones donde coincidan los nombres tienen LCR = t
    LCR$bold_vs_midori[ which(LCR[, midori_rank] == LCR[, bold_rank]) ] <- t
}

# bold vs coarbitrator ----

LCR$bold_vs_coarbitrator <- 1
# recorrer hacia el 7 en loop
for(t in 1:7){
  # Generar nombre de la columna con taxonomia por morfologia
  coarbitrator_rank <- paste("coarbitrator", TL[t], sep="_")
  # Generar nombre de la columna con taxonomica por molecular
  bold_rank <- paste("bold", TL[t], sep="_")
  # Aquellos renglones donde coincidan los nombres tienen LCR = t
  LCR$bold_vs_coarbitrator[ which(LCR[, midori_rank] == LCR[, bold_rank]) ] <- t
}


# getback the Last Common Rank

bold_LCR <-NULL
for (i in 1:nrow(bold)) {
    rl <- LCR[,'bold_SL'][i] + 1
    bold_LCR[[i]] <- list(rank=names(bold)[rl], linage=bold[i,rl]) }
bold_LCR <- do.call(rbind, bold_LCR) 
bold_LCR <- data.frame(rank = do.call(rbind, bold_LCR[,1]), lineage = do.call(rbind, bold_LCR[,2]))
bold_LCR$rank <- sapply(strsplit(as.vector(bold_LCR$rank), "_"), "[[", 2)

# ES NECESARIO VERIFICAR ESTA PARTE,
midori_LCR <-NULL
for (i in 1:nrow(midori)) {
    rl <- LCR[,'midori_SL'][i] + 1
    midori_LCR[[i]] <- list(rank=names(midori)[rl], linage=midori[i,rl]) }
midori_LCR <- do.call(rbind, midori_LCR)
midori_LCR <- data.frame(rank = do.call(rbind, midori_LCR[,1]), lineage = do.call(rbind, midori_LCR[,2]))
midori_LCR$rank <- sapply(strsplit(as.vector(midori_LCR$rank), "_"), "[[", 2)

coarbitrator_LCR <-NULL
for (i in 1:nrow(coarbitrator)) {
  rl <- LCR[,'bold_vs_coarbitrator'][i] + 1
  coarbitrator_LCR[[i]] <- list(rank=names(coarbitrator)[rl], linage=coarbitrator[i,rl]) }
coarbitrator_LCR <- do.call(rbind, coarbitrator_LCR)
coarbitrator_LCR <- data.frame(rank = do.call(rbind, coarbitrator_LCR[,1]), lineage = do.call(rbind, coarbitrator_LCR[,2]))
coarbitrator_LCR$rank <- sapply(strsplit(as.vector(coarbitrator_LCR$rank), "_"), "[[", 2)

# and save data imput
data <- as_tibble(data.frame( ASV = midori$ASV,
                    #midori = do.call(rbind, midori_[,1]), 
                    #bold = do.call(rbind, bold_[,1]),
                    lineage_x = do.call(rbind, bold_[,2]),
                    lineage_y = do.call(rbind, midori_[,2]), 
                    lineage_z = do.call(rbind, coarbitrator_[,2]),
                    # Anadimos el valor de cambio en las asignaciones:
                    rank_x = bold[,9],
                    rank_y = midori[,9], 
                    rank_z = coarbitrator[,9],
                    diff_x_y =  bold[,9] - midori[,9], 
                    diff_x_z =  bold[,9] - coarbitrator[,9],
                    LCR_x_y = LCR$bold_vs_midori,
                    LCR_x_z = LCR$bold_vs_coarbitrator,
                    LCR_x = bold_LCR$lineage,
                    LCR_y = midori_LCR$lineage,
                    LCR_z = coarbitrator_LCR$lineage,
                    stringsAsFactors =  FALSE)
                    )  # NOT ABS VALUE IN DIFF

data$diff_x_y_z <- data$diff_x_y - data$diff_x_z

# Calculate distance:

# Computes pairwise string distances between elements
library('stringdist')

# We  use is the Levenshtein Distance which is simply the 
# single character edit distance between two strings by 
# a combination of insertions, deletions, or substitutions

dist_x_y <- stringdist(data$lineage_x, data$lineage_y, method = 'lv') # inf values reported
dist_x_z <- stringdist(data$lineage_x, data$lineage_z, method = 'lv') # es necesario arreglar la nomenclatura de coarbitrator

# esto cambia

filter(dist_x_y, dist > 0) %>%
  select(midori, bold, diff) %>%
  melt(id.vars = "diff", variable.name = "DB", value.name = "Rank") %>% 
  as_tibble() %>%   # mutate(diff = abs(diff))
  mutate(Rank = factor(Rank, levels = TL[-1]) ) -> y
# ahora comprobamos que los niveles de cambio zero sean consistentes en su asignacion:

disvis <- dist[dist$diff == 0 & dist$dist > 0,]
disvis <- dist[ dist$dist > 0,]

filter(dist, dist > 0) %>%
  select(LCR, LCR_x, LCR_y, dist) %>%
  melt(id.vars = c("LCR", "dist"), variable.name = "LCR_xy", value.name = "lineage") %>%
  as_tibble() -> dist_xy

p <- ggplot(dist_xy, aes(x=dist, fill=lineage)) +
  #geom_histogram(position="dodge2", bins = 40) +
  geom_density(alpha = 0.6) +
  facet_wrap (~ lineage , scales = 'free_y') +
  theme(legend.position="top")
p

p <- ggplot(dist_xy, aes(x=dist)) +
  geom_density(alpha = 0.6) + facet_wrap (~ LCR_xy, scales = 'free_y')


plot(density(dist$dist), main = 'Comparing DB assignations by String Distance [LCR]')
lines(stats::lowess(dist$dist))

d <- stringdistmatrix(disvis$LCR_x, disvis$LCR_y, useNames="none", method = 'lv') 
library(superheat)

superheat(d, row.dendrogram = TRUE, col.dendrogram = TRUE)

quit(save = 'no')

data$weight <- 'None'
data[data$diff > 0 , 'weight'] <- 'midori'
data[data$diff < 0 , 'weight'] <- 'bold'

filter(data, weight != 'None') %>%
  select(midori, bold, diff) %>%
  melt(id.vars = "diff", variable.name = "DB", value.name = "Rank") %>% 
  as_tibble() %>%   # mutate(diff = abs(diff))
  mutate(Rank = factor(Rank, levels = TL[-1]) ) -> y


# peces entre las bases de datos (bold)
ntaxa <-  function(r) {
  r <- data.frame(table(r))  
  r <- r[order(r$Freq, decreasing = TRUE), ]
  # r <- r[order(r$r, decreasing = FALSE), ]
  return(r) }


ntaxa(bold[,4]) # Actinopterygii, 691
ntaxa(midori[,4]) # Actinopteri, 582
ntaxa(coarbitrator[,4]) # Actinopteri, 388

# arreglar la taxonomia primero
acti_bold <- midori[bold[,4] == 'Actinopterygii', ] # 13056
acti_midori <- midori[midori[,4] == 'Actinopteri', ] # 13014
acti_coarbitrator <- coarbitrator[coarbitrator[,4] == 'Actinopteri', ] # 13023




# Barplot
# use SL_diff or SL_tidy
# Midori and bold taxonomy resolution
# Won (+) and lost (−) levels (shoul)

SL <- select(data, ASV, diff_x_y, diff_x_z) 
SL$weight_x_y <- 'shared_x_y'
SL[SL$diff_x_y > 0 , 'weight_x_y'] <- 'midori'
SL[SL$diff_x_y < 0 , 'weight_x_y'] <- 'bold'

SL$weight_x_z <- 'shared_x_z'
SL[SL$diff_x_z > 0 , 'weight_x_z'] <- 'coarbitrator'
SL[SL$diff_x_z < 0 , 'weight_x_z'] <- 'bold'

SL[,'diff_x_y'] <- abs(SL[,'diff_x_y'])
SL[,'diff_x_z'] <- abs(SL[,'diff_x_z'])

table(SL$weight_x_z)

SL_melt <- melt(SL, 
                id.vars= c("ASV","diff_x_y", "diff_x_z"),
                value.name = "DB",
                variable.name = 'weight')

SL_melt <- melt(SL_melt)

SL_diff <- aggregate(SL_melt[,'value'], by=list(SL_melt[,'DB']), FUN = table)


SL_diff <- data.frame(level = SL_diff[,1], SL_diff[,2])
SL_diff <- melt(SL_diff)
# barplot


ggplot(SL_melt, aes(value, variable, fill = diff)) + 
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  geom_text(aes(label=ASVs), size=3, vjust=-0.2) +
  scale_fill_brewer(palette = "Set1") #+
  #labs(x="Levels changed", y="ASVs (n)", 
       title="Midori - Bold taxonomy Difference", 
       subtitle=paste0("Midori and Bold levels\nTotal ASVs = ", nrow(midori), "\n",
                       "ASVs without levels change across data bases: ",no_change),
       caption=paste0("Number of ASVs with levels change across databases: ","\n",
                      "Midori: ", midori_lc, " and Bold: ", bold_lc)) +
  theme()
# ALLUVIAL

#filter(data, weight != 'shared') %>%
data %>%
  select(rank_x, rank_y, rank_z) %>%
  mutate(diff = abs(diff), weight = as.factor(weight),
         midori = factor(midori, levels = TL[-1]), 
         bold = factor(bold, levels = TL[-1])) -> alluvia

# library('ggalluvial')
ggplot(data = alluvia,
       aes(axis1 = midori, axis2 = bold)) +
  geom_alluvium(aes(fill = weight)) + 
  geom_stratum() + 
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_text(stat = "stratum", label.strata = TRUE) +
  # theme(legend.position = "bottom") +
  theme_minimal() 

library(alluvial)

alluvial(select(data, lineage_x, lineage_y, lineage_z), freq=data$diff_x_y_z,
         col = ifelse(data$diff_x_y_z > 0, "orange", "grey"),
         border = ifelse(data$diff_x_y_z > 0, "orange", "grey"),
         hide = data$diff_x_y_z == 0,
         cex = 0.7)

head(as.data.frame(Titanic, stringsAsFactors = FALSE))

#
# Parallel Coordinates
library(GGally)

paralle <- select(data, rank_x, rank_y, rank_z, diff_x_y, diff_x_z)

# sampleTree = hclust(dist(paralle[,1:3]), method = "average");
# plot(sampleTree, main = "", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

cbind(data.frame(table(paralle$diff_x_y)), data.frame(table(paralle$diff_x_z)))

paralle$weight_x_y <- 'shared_x_y'
paralle[paralle$diff_x_y > 0 , 'weight_x_y'] <- 'midori'
paralle[paralle$diff_x_y < 0 , 'weight_x_y'] <- 'bold'

paralle$weight_x_z <- 'shared_x_z'
paralle[paralle$diff_x_z > 0 , 'weight_x_z'] <- 'coarbitrator'
paralle[paralle$diff_x_z < 0 , 'weight_x_z'] <- 'bold'

ggparcoord(data, columns = 4:8, groupColumn = "sex")

ggparcoord(paralle, columns = 1:3, groupColumn = "weight_x_y")



# Lineplot. Distribucion de las asignaciones entre bases de datos
L_tables <- data.frame(Level=-5:5, 
                       x_y=data.frame(table(data$diff_x_y))[,'Freq'], 
                       x_z=data.frame(table(data$diff_x_z))[, 'Freq'])

L_tidy <- melt(L_tables, id.vars = 'Level',
               variable.name = 'diff',
               value.name = 'n')
# wrap by standar (BOLD) and reference (midori and GBank)
L_tidy$Ref <- 'Standard(BOLD)'
L_tidy[L_tidy$Level > 0, 'Ref'] <- 'Reference'

# 
L_tidy$database <- ''
L_tidy[L_tidy$diff == 'x_y', 'database'] <- 'Midori'
L_tidy[L_tidy$diff == 'x_z', 'database'] <- 'GBank'    
       
L_tidy$Level <- abs(L_tidy$Level)


L_tidy$Level <- factor(L_tidy$Level, levels = 0:5, ordered = T)



ggplot(subset(L_tidy, Level > 0), aes(x=Level, y=n, group=database, fill=database, color=database)) + 
  geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") +
  geom_text(data = subset(L_tidy, Level > 0), aes(label=n), size=4, vjust=1, color = 'black') +
  # geom_text(aes(label=diff), size=3, vjust=-0.2, color = 'black') +
  scale_color_brewer(palette = "Set1") + 
  #theme(panel.grid = element_blank(),
  #      panel.background = element_blank()) +
  labs(y="ASVs", title="Taxonomy resolution Distribution",
       subtitle=paste0("Total ASVs in (run14) = ", nrow(midori), "\n",
                       "Undetermined Bold = ", n_undetermined_bold, "\n",
                       "Undetermined Midori = ", n_undetermined_midori, "\n",
                       "Undetermined GBank = ", n_undetermined_coarbitrator),
       caption = 'Reference panel show how best the database was against the standard \nwhereas the Standard(BOLD) panel shows how best the BOLD db was against either, Midori or GBank db') +
  facet_grid(~Ref)

