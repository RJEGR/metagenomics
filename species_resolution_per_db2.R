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


SL_diff_ <- data.frame(table(data$diff), stringsAsFactors = FALSE)

# apply(data[, c('midori', 'bold')], 2, table)
y <- apply(data[data$diff != 0, c('midori', 'bold')], 2, table)
# Order taxa by TL
y <- y[match(TL[-1], rownames(y)),]

# SL_data <- as.data.frame(data)
# SL_data_agg <- aggregate(SL_data[, c('midori', 'bold')], by=list(SL_data[,'diff']), FUN = table)

data$weight <- 'None'
data[data$diff > 0 , 'weight'] <- 'midori'
data[data$diff < 0 , 'weight'] <- 'bold'

filter(data, weight != 'None') %>%
  select(midori, bold, diff) %>%
  melt(id.vars = "diff", variable.name = "DB", value.name = "Rank") %>% 
  as_tibble() %>%   # mutate(diff = abs(diff))
  mutate(Rank = factor(Rank, levels = TL[-1]) ) -> y

#filter(y, diff > 0 & DB == 'midori') %>% as.data.frame() -> test
#filter(y, diff < 0 & DB == 'bold') %>% as.data.frame() -> test


#y_agg <- aggregate(test[, c('Rank', 'DB')], by=list(test[,'diff']), FUN = table)
#y_agg <- t(data.frame(diff = y_agg[,1], y_agg[,2]))[-1,]
#y_agg <- y_agg[match(TL[-1], rownames(y_agg)),]




# Somo alluvials


filter(data, weight != 'None') %>%
      select(midori, bold, diff, weight) %>%
      mutate(diff = abs(diff), weight = as.factor(weight),
              midori = factor(midori, levels = TL[-1]), 
              bold = factor(bold, levels = TL[-1])) -> alluvia

library('ggalluvial')

ggplot(data = alluvia,
       aes(axis1 = midori, axis2 = bold)) +
  geom_alluvium(aes(fill = weight)) + 
  geom_stratum() + 
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_text(stat = "stratum", label.strata = TRUE) +
  # theme(legend.position = "bottom") +
  theme_minimal() 


# Computes pairwise string distances between elements
library('stringdist')

# We  use is the Levenshtein Distance which is simply the 
# single character edit distance between two strings by 
# a combination of insertions, deletions, or substitutions

dist <- as_tibble(data.frame(data, dist = stringdist(data$lineage_x, data$lineage_y, method = 'lv'))) # inf values reported


# ahora comprobamos que los niveles de cambio zero sean consistentes en su asignacion:

disvis <- dist[dist$diff != 0 ,]


plot(density(disvis$dist), main = 'Comparing DB assignations by String Distance')
lines(stats::lowess(disvis$dist))


d <- stringdistmatrix(disvis$lineage_x, disvis$lineage_y, useNames="none", method = 'lv') 
library(superheat)

superheat(d, row.dendrogram = TRUE, col.dendrogram = TRUE)

#d <- stringdistmatrix(disvis$lineage_x, disvis$lineage_y, useNames="strings")
#h <- stats::hclust(lineage_matrix)
#plot(d)

quit(save = 'no')

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

# Barplot

ggplot(SL_diff_, aes(Var1, Freq, fill)) + 
    geom_col(width = 0.4, fill="thistle4") +
    geom_text(aes(label=Freq), size=3, vjust=-0.2) +
    labs(x="Levels changed", y="ASVs (n)", 
        title="Midori - Bold taxonomy Difference", 
        subtitle=paste0("Midori (+) and Bold (-) levels\nTotal ASVs = ", nrow(midori))) +
    theme()

# use SL_diff or SL_tidy
# Midori and bold taxonomy resolution
# Won (+) and lost (−) levels (shoul)

SL <- data
# SL$Reward <- ifelse(SL$diff > 0, "Midori",
#                     ifelse(SL$diff < 0 ,"Bold", "equal"))

SL[,'diff'] <- abs(SL[,'diff'])

SL <- as.data.frame(SL)

SL_diff <- aggregate(SL[,'weight'], by=list(SL[,'diff']), FUN = table)

SL_diff <- data.frame(level = SL_diff[,1], SL_diff[,2])

SL_tidy <- melt(SL_diff, id.vars="level",
                         variable.name = "DB", 
                         value.name = "ASVs")

SL_tidy$level <- abs(SL_tidy$level)
# SL_tidy$Rank <- factor(SL_tidy$rank, levels= TL[2:7], labels= TL[2:7], ordered = T)

datavis <- SL_tidy[SL_tidy$DB != 'None' & SL_tidy$ASVs !=0,]
no_change <- SL_tidy[SL_tidy$DB == 'None' & SL_tidy$ASVs !=0, 'ASVs']
bold_lc <- sum(SL_tidy[SL_tidy$DB != 'None' & SL_tidy$DB == 'bold', 'ASVs'])
midori_lc <- sum(SL_tidy[SL_tidy$DB != 'None' & SL_tidy$DB == 'midori', 'ASVs'])
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

