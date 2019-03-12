#!/usr/bin/env Rscript

# Load libraries ----
.cran_packages <- c("stringr", "reshape2") 
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

bold.file = "run012_relax_ASVs.BOLD_public_species.wang.taxonomy_99"
midori.file  = "run012_relax_ASVs.midori_unique_DB_0.wang.taxonomy_99"



# To load file:  wang.taxonomy
source(file = "~/Documents/GitHub/metagenomics/converted2worms/cigom_load_dada2_tables.R")

## Load assignation result ----

bold <- read.wangtax(bold.file)
midori <- read.wangtax(midori.file)

if (identical(names(bold), names(midori))) { TL <- names(bold)[-c(1,ncol(bold))]
    } # else {TL <- rank}

# Add column for converted level (CL)

# plot(density(rdp_df$SL))

n_undetermined_bold <- nrow(bold[bold$Phylum == 'Undetermined',])
n_undetermined_midori <- nrow(midori[midori$Phylum == 'Undetermined',])

# Plots ----


x <- data.frame(Rank=TL[-1],
            midori=data.frame(table(midori$SL))[,2],
            bold=data.frame(table(bold$SL))[,2])

midori_ <- NULL
for (i in 1:nrow(midori)) {
    rl <- midori$SL[i] + 1
    midori_[[i]] <- list(rank=names(x)[rl], linage=x[i,rl])
    
}

bold_ <-NULL
for (i in 1:nrow(bold)) {
    rl <- bold$SL[i] + 1
    bold_[[i]] <- list(rank=names(x)[rl], linage=x[i,rl])
}



midori_ <- do.call(rbind, midori_)
bold_ <- do.call(rbind, bold_)

data <- data.frame( ASV = midori$ASV,
                    midori = do.call(rbind, midori_[,1]), 
                    lineage_x = do.call(rbind, midori_[,2]),
                    bold = do.call(rbind, bold_[,1]),
                    lineage_y = do.call(rbind, bold_[,2]))



data %>% 
  distinct(lineage_x, lineage_y) ->  lineages 

data[,c(2,3)] %>% group_by(midori) %>%
  summarise(n = table(lineage_x)) -> test


data_tidy <- melt(data, id.vars = c("ASV", "lineage_x", "lineage_y"),
                    variable.name = "DB",
                    value.name = "Rank")

# or

data_tidy <- melt(data[-1], id.vars=c("midori", "bold"))

# aggregate(data_tidy[,'DB'], by=list(data_tidy[,'Rank']), FUN = table)

data_agg <- aggregate(data_tidy[,c('lineage_x', 'lineage_y')], by=list(data_tidy[,'Rank']), FUN = table)
data_agg <- data.frame(level = SL_diff[,1], SL_diff[,2])
# bastante complicado elaborar el formato del sanky pplot

L_tidy$Rank <- factor(L_tidy$Rank, levels = TL[2:7], ordered = T)

testvis <- data_tidy[sample(nrow(data_tidy), 100), ]

ggplot(data = x,
       aes(axis1 = midori, axis2 = bold)) +
  geom_alluvium(aes(fill = Rank)) + #aes(fill = DB)
  geom_stratum() + 
  geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("")



# BOLD_public_species table(bold$bold_SL)
#   2    3    4    5    6    7 
# 8567  176  198  322  352 4257 
# Removing 7937 Undetermined Phylum from the dataset

# midori
#   2    3    4    5    6    7 
# 8666   96  214  308  396 4192 
# # Removing 8000 Undetermined Phylum from the dataset
# make function to read an process by wangtax

# Lineplot. Distribucion de las asignaciones entre bases de datos
L_tables <- data.frame(Rank=TL[2:7], Level=2:7, 
                       bold=data.frame(table(bold$SL))[,'Freq'], 
                       midori=data.frame(table(midori$SL))[, 'Freq'])

L_tables$diff <- abs(L_tables$bold - L_tables$midori)

# tidy data (reshape2)
L_tidy <- melt(L_tables, id.vars=c("Rank", "Level", "diff"),
                variable.name = "DB",
                value.name = "ASVs")
L_tidy$Rank <- factor(L_tidy$Rank, levels = TL[2:7], ordered = T)


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
# SL_tidy$Rank <- factor(SL_tidy$rank, levels= TL[2:7], labels= TL[2:7], ordered = T)

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

# Basicamente se necesita para cada ASV la taxonomia (una columna por nivel) con DB_1 y con DB_2

names(midori)[-1] <- paste("midori", names(midori)[-1], sep="_")
names(bold)[-1] <- paste("bold", names(bold)[-1], sep="_")

input_midori_bold <- cbind(midori, bold[,-1])

# aniadir una columna, rango 1
input_midori_bold$midori_vs_bold <- 1
# recorrer hacia el 7 en loop
for(t in 1:7){
    # Generar nombre de la columna con taxonomia por morfologia
    midori_rank <- paste("midori", TL[t], sep="_")
    # Generar nombre de la columna con taxonomica por molecular
    bold_rank <- paste("bold", TL[t], sep="_")
    # Aquellos renglones donde coincidan los nombres tienen LCR = t
    input_midori_bold$midori_vs_bold[ which(input_midori_bold[,midori_rank] == input_midori_bold[,bold_rank]) ] <- t
}

# t(input_midori_bold[11,])
# paste("midori", TL[-1], sep="_")



input_midori_bold_m <- melt(input_midori_bold, id.vars=c("ASV", "midori_SL", "bold_SL", "midori_vs_bold"),
                variable.name = "Rank",
                value.name = "linage")

input_midori_bold_m$DB <- sapply(strsplit(as.vector(input_midori_bold_m$Rank), "_"), "[[", 1)
input_midori_bold_m$Rank <- sapply(strsplit(as.vector(input_midori_bold_m$Rank), "_"), "[[", 2)


sankey_data <- input_midori_bold_m[sample(nrow(input_midori_bold_m), 100), ] # to test

sankey_data %>% group_by(Rank, linage, DB) %>%
  summarise(n = sum(midori_vs_bold)) -> sankey_vis

ggplot(data = sankey_vis,
       aes(axis1 = Rank, axis2 = linage,
           y = n)) +
  #scale_x_discrete(limits = c("Class", "Sex", "Age"), expand = c(.1, .05)) +
 # xlab("Demographic") +
  geom_alluvium(aes(fill = DB)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("")






# Percentage of correlacion entre midori y bold
midori_bold_corr <- 100-round((table(input_midori_bold$midori_vs_bold)/nrow(input_midori_bold))*100, 2)
midori_bold_corr

# match_rank / n_asvs * 100

cor <- input_midori_bold[, c("midori_SL", "bold_SL", "midori_vs_bold")]
midori_vs_bold <- input_midori_bold[input_midori_bold$midori_vs_bold == 7, c('bold_Species', 'midori_Species')]


ggplot(table_m, 
       aes(x=item_morfo, y=rank_txt, shape=resolution_type, size=resolution_type, colour=resolution_type)) +
    geom_point()

table_m[which(table_m$Query_ID %in% "ConcensoIctio25"), ]


#          Query_ID                              item_morfo morfo_Order resolution_type rank_txt            solo_morfo
# 1  ConcensoIctio25 ConcensoIctio25\nVinciguerria poweriae      Stomii       Molecular    Genus Vinciguerria poweriae
# 34 ConcensoIctio25 ConcensoIctio25\nVinciguerria poweriae      Stomii   Morphological  Species Vinciguerria poweriae
# 67 ConcensoIctio25 ConcensoIctio25\nVinciguerria poweriae      Stomii   Matching rank    Genus Vinciguerria poweriae

# Plot
poster_plot <- ggplot(table_m, 
       aes(x=item_morfo, y=rank_txt, shape=resolution_type, size=resolution_type, colour=resolution_type)) +
    geom_point() +
    facet_grid(morfo_Order ~., scales="free", space="free") +
    labs(title="", x="", y="") +
    scale_shape_manual(name="Identification approach", values=c(16, 1, 2)) +
    scale_colour_manual(name="Identification approach", values=color_asign) +
    scale_size_manual(values=c(6,7.2,7.5)) +
    guides(size=F, colour = guide_legend(override.aes = list(size=6))) +
    coord_flip() +
    theme_light() +
    theme(axis.text.x = element_text(angle=0, size=12, hjust=0.5, vjust=0.5),
          axis.text.y = element_text(size=12),
          strip.text.y = element_text(angle=0, size=12),
          #strip.background = element_rect(colour="grey40", fill=NA),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12, face="bold"),
          legend.key = element_rect(fill="white"),
          legend.position = "top")
poster_plot  


#alluvial
# devtools::install_github("r-lib/rlang", build_vignettes = TRUE)
detach("package:tidyverse", unload=TRUE)
library(ggalluvial)
titanic_wide <- data.frame(Titanic)
head(titanic_wide)

ggplot(data = titanic_wide,
       aes(axis1 = Class, axis2 = Sex, axis3 = Age,
           y = Freq)) +
  scale_x_discrete(limits = c("Class", "Sex", "Age"), expand = c(.1, .05)) +
  xlab("Demographic") +
  geom_alluvium(aes(fill = Survived)) +
  geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE) +
  theme_minimal() +
  ggtitle("passengers on the maiden voyage of the Titanic",
          "stratified by demographics and survival")


head(sp <- data.frame(table(bold[!is.na(bold$bold_Species) ,8])))

sp[c(grep(' sp[.]', sp[,1])),]

input_midori_bold_m[input_midori_bold_m$Rank == 'Species' & input_midori_bold_m$]