#!/usr/bin/env Rscript

# Load libraries ----
.cran_packages <- c("stringr", "tidyverse", "reshape2", "ggplot2") 


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

bold.file = "run012_relax_ASVs.BOLD_public_species.wang.taxonomy_99"
midori.file  = "run012_relax_ASVs.midori_unique_DB_0.wang.taxonomy_99"



# To load file:  wang.taxonomy
source(file = "~/Documents/GitHub/metagenomics/converted2worms/cigom_load_dada2_tables.R")

## Load assignation result ----

"read.wangtax" <- function(x) {
    
    rdp_df <- NULL
# wang.taxonomy, threshold = 99
    rdp_df <- load.mothurdp_boots(x, 80)  #output is data.frame 'rdp_df'

# Remove unclassified and replace underscores by spaces
    rdp_df <- mutate_all(data.frame(rdp_df), funs(str_replace_all(., c(".+_unclassified"=NA, "_"=" "))))

    rdp_df$ASV <- str_replace(rdp_df$ASV, " ", "_")

    max.rank <-  ncol(rdp_df) -1 

# Tag starting level (SL) (aka calculate resolution)
    rdp_df$SL <- max.rank

    for(t in max.rank:2){
       rr <- t-1  #rank real 
       rc <- t+1  #column-index (rank) para corregir
            if(t == max.rank){
               rdp_df$SL[is.na(rdp_df[,rc])] <- rr
                } else {
                rdp_df$SL[is.na(rdp_df[,rc]) & rdp_df$SL <= t ] <- rr
        }
    }

return(rdp_df) 
}


bold <- read.wangtax(bold.file)
midori <- read.wangtax(midori.file)

if (identical(names(bold), names(midori))) { TL <- names(bold) 
    } # else {TL <- rank}

# Add column for converted level (CL)

# plot(density(rdp_df$SL))

n_undetermined_bold <- nrow(bold[bold$Phylum == 'Undetermined',])
n_undetermined_midori <- nrow(midori[midori$Phylum == 'Undetermined',])

# Plots ----

# BOLD_public_species
#   2    3    4    5    6    7 
# 8567  176  198  322  352 4257 
# Removing 7937 Undetermined Phylum from the dataset

# midori
#   2    3    4    5    6    7 
# 8666   96  214  308  396 4192 
# # Removing 8000 Undetermined Phylum from the dataset
# make function to read an process by wangtax

# Lineplot. Distribucion de las asignaciones entre bases de datos
L_tables <- data.frame(Rank=TL[3:8], Level=2:7, 
                       bold=data.frame(table(bold$SL))[,'Freq'], 
                       midori=data.frame(table(midori$SL))[, 'Freq'])

L_tables$diff <- abs(L_tables$bold - L_tables$midori)

# tidy data (reshape2)
L_tidy <- melt(L_tables, id.vars=c("Rank", "Level", "diff"),
                variable.name = "DB",
                value.name = "ASVs")
L_tidy$Rank <- factor(L_tidy$Rank, levels = TL[3:8], ordered = T)


# or with tidy
# L_tidy <- gather(L_tables, key="DB", value = "ASVs", -c("Rank", "Level"))
# L_tables_g$Stage <- factor(L_tables_g$DB, levels=c("bold", "midori"), labels=c("Initial", "Converted"), ordered = T)


ggplot(L_tidy, aes(x=Rank, y=ASVs, group=DB, fill=DB, color=DB)) + 
    geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") +
    geom_text(data = subset(L_tidy, DB == 'midori'), aes(label=diff), size=3, vjust=0.2, color = 'black') +
    # geom_text(aes(label=diff), size=3, vjust=-0.2, color = 'black') +
    scale_color_brewer(palette = "Set1") + 
    #theme(panel.grid = element_blank(),
    #      panel.background = element_blank()) +
    labs(x="", title="Midori - Bold Taxonomy resolution Distribution",
               subtitle=paste0("Total ASVs = ", nrow(midori), "\n",
                               "Undetermined midori = ", n_undetermined_midori, "\n",
                               "Undetermined bold = ", n_undetermined_bold))

# n_undetermined_midori

# Starting level (SL) (aka calculate resolution)
SL <- data.frame(bold = bold$SL, 
                 midori = midori$SL, 
                 diff =  midori$SL - bold$SL) #not abs() values

rownames(SL) <- rdp_df$ASV
# Midori and bold taxonomy resolution

# Won (+) and lost (−) levels (shoul)

SL$Reward <- ifelse(SL$diff > 0, "Midori",
                      ifelse(SL$diff < 0 ,"Bold", "equal"))

table(SL$Reward)

SL_diff_ <- data.frame(table(SL$diff), stringsAsFactors = FALSE)

# Barplot

ggplot(SL_diff_, aes(Var1, Freq, fill)) + 
    geom_col(width = 0.4, fill="thistle4") +
    geom_text(aes(label=Freq), size=3, vjust=-0.2) +
    labs(x="Levels changed", y="ASVs (n)", 
        title="Midori - Bold taxonomy Difference", 
        subtitle=paste0("Midori (+) and Bold (-) levels\nTotal ASVs = ", nrow(midori))) +
    theme()

# use SL_diff or SL_tidy

SL_diff <- aggregate(SL[,'Reward'], by=list(SL[,'diff']), FUN = table)

SL_diff <- data.frame(level = SL_diff[,1], SL_diff[,2])

SL_tidy <- melt(SL_diff, id.vars="level",
                         variable.name = "DB", 
                         value.name = "ASVs")

SL_tidy$level <- abs(SL_tidy$level)
# SL_tidy$Rank <- factor(SL_tidy$rank, levels= TL[3:8], labels= TL[3:8], ordered = T)

datavis = SL_tidy[SL_tidy$DB != 'equal' & SL_tidy$ASVs !=0,]
no_change = SL_tidy[SL_tidy$DB == 'equal' & SL_tidy$ASVs !=0, 'ASVs']
bold_lc <- sum(SL_tidy[SL_tidy$DB != 'equal' & SL_tidy$DB == 'Bold', 'ASVs'])
midori_lc <- sum(SL_tidy[SL_tidy$DB != 'equal' & SL_tidy$DB == 'Midori', 'ASVs'])
# barplot

ggplot(datavis, aes(level, ASVs, fill = DB)) + 
    geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
    geom_text(aes(label=ASVs), size=3, vjust=-0.2) +
    scale_fill_brewer(palette = "Set1") +
    labs(x="Levels changed", y="ASVs (n)", 
        title="Midori - Bold taxonomy Difference", 
        subtitle=paste0("Midori and Bold levels\nTotal ASVs = ", nrow(midori), "\n",
                        "ASVs without levels change across data bases: ",no_change),
        caption=paste0("Number of ASVs with levels change across databases: ","\n",
                        "Midori: ", midori_lc, " and Bold: ", bold_lc)) +
    theme()

quit(save = 'no')





 Bold  equal Midori 
   577  12795    500 

   Es decir, hay una diferencia de 500/577 mejor asignaciones entre midori y bold, mientras que 12,795 ASVs
   alcanzaron la misma resolucion usando ambas bases de referencia

        "bold"   "midori" "diff"   "Reward"
ASV_13745    7      2   -5   Bold
ASV_13778    7      6   -1   Bold
ASV_13852    2      2    0  equal
ASV_13853    3      6    3 Midori
ASV_13854    2      7    5 Midori

