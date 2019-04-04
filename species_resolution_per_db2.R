#!/usr/bin/env Rscript

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

bold.file = "run012_relax_ASVs.BOLD_public_species.wang.taxonomy_99"
midori.file  = "run012_relax_ASVs.midori_unique_DB_0.wang.taxonomy_99"

# To load file:  wang.taxonomy
source(file = "~/Documents/GitHub/metagenomics/converted2worms/cigom_load_dada2_tables.R")

## Load assignation result ----

bold <- read.wangtax(bold.file)
midori <- read.wangtax(midori.file)

if (identical(names(bold), names(midori))) { TL <- names(bold)[-c(1,ncol(bold))]
    } 

# Contamos el numero de indeterminados por base de datos a nivel filo
n_undetermined <- function(x) ( nrow(x[x[,3] == 'Undetermined',]) )
n_undetermined_bold <- n_undetermined(bold)
n_undetermined_midori <- n_undetermined(midori)

# Recuperamos el tag (asignacion) de cada una de las posiciones en SL

x <- data.frame(Rank=TL[-1],
            midori=data.frame(table(midori$SL))[,2],
            bold=data.frame(table(bold$SL))[,2])

# Back the SL position per database
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

# compare the last parent between db

names(midori)[-1] <- paste("midori", names(midori)[-1], sep="_")
names(bold)[-1] <- paste("bold", names(bold)[-1], sep="_")

LCR <- data.frame(midori, bold[,-1], diff =  midori[,9] - bold[,9], stringsAsFactors =  FALSE)
# aniadir una columna, rango 1
LCR$midori_vs_bold <- 1
# recorrer hacia el 7 en loop
for(t in 1:7){
    # Generar nombre de la columna con taxonomia por morfologia
    midori_rank <- paste("midori", TL[t], sep="_")
    # Generar nombre de la columna con taxonomica por molecular
    bold_rank <- paste("bold", TL[t], sep="_")
    # Aquellos renglones donde coincidan los nombres tienen LCR = t
    LCR$midori_vs_bold[ which(LCR[, midori_rank] == LCR[, bold_rank]) ] <- t
}
#

#input_midori_bold_m <- melt(input_midori_bold, id.vars=c("ASV", "midori_SL", "bold_SL", "midori_vs_bold"),
#                variable.name = "Rank",
#                value.name = "linage")

#input_midori_bold_m$DB <- sapply(strsplit(as.vector(input_midori_bold_m$Rank), "_"), "[[", 1)
# input_midori_bold_m$Rank <- sapply(strsplit(as.vector(input_midori_bold_m$Rank), "_"), "[[", 2)



data <- as_tibble(data.frame( ASV = midori$ASV,
                    midori = do.call(rbind, midori_[,1]), 
                    bold = do.call(rbind, bold_[,1]),
                    lineage_x = do.call(rbind, midori_[,2]),
                    lineage_y = do.call(rbind, bold_[,2]), 
                    # Anadimos el valor de cambio en las asignaciones:
                    rank_x = midori[,9],
                    rank_y = bold[,9], 
                    diff =  midori[,9] - bold[,9], 
                    LCR = LCR$midori_vs_bold,
                    stringsAsFactors =  FALSE)
                    )  # NOT ABS VALUE IN DIFF

quit(save = 'no')