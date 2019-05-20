#!/usr/bin/env Rscript

rm(list=ls()); 

# Functions ----



mrank <- function(x) {
  x_ <- NULL
  for (i in 1:nrow(x)) {
  rl <- x$SL[i] + 1
  # x_[[i]] <- list(rank=names(x)[rl], linage=x[i,rl]) }
  x_[[i]] <- names(x)[rl]
  }
}

source(file = "~/Documents/GitHub/metagenomics/readtx.R")

library(ggplot2)
library(reshape2)
# Set filenames ----
path_BOLD <- '/Users/cigom/metagenomics/COI/species_resolution_per_db'
bold_all <- 'run014_t2_ASVs.ALL.wang.taxonomy'
bold_sp <- 'run014_t2_ASVs.BOLD_public_species.wang.taxonomy'


# Load files ----
full.taxa.obj <- read_rdp(paste0(path_BOLD,'/',bold_all))
sp.taxa.obj <- read_rdp(paste0(path_BOLD,'/',bold_sp))

TL <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

colnames(full.taxa.obj) <- c(TL, 'SL')
colnames(sp.taxa.obj) <- c(TL, 'SL')



# Level Resolution vis ----

SL <- rbind(data.frame(SL = full.taxa.obj$SL,
                 DataBase = 'full'),
            data.frame(SL = sp.taxa.obj$SL,
                       DataBase = 'SP'))


# SL <- SL[SL$SL > 0 ,]
SL$SL <- factor(SL$SL, levels = names(table(SL$SL)))

# and
SL_agg <- aggregate(SL[,'DataBase'], by=list(SL[,'SL']), FUN = table)
SL_agg <- data.frame(level = SL_agg[,1], SL_agg[,2])
SL_agg$Rank <- TL[-1]
SL_aggM <- melt(SL_agg)
SL_aggM$Rank <- factor(SL_aggM$Rank, levels = TL[-1])
SL_aggM$variable <- factor(SL_aggM$variable, levels = c("full","SP"))

ggplot(SL_aggM, aes(x=Rank, y=value, group=variable, fill=variable, color=variable)) + 
  geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") +
  geom_text(data = SL_aggM, aes(label=value), size=4, vjust=1, color = 'black') +
  scale_color_brewer(palette = "Set1") +
  labs(title="BOLD species & BOLD full db taxonomy Difference",
       x = 'Taxonomy level',
       y = '# ASVs')
       #subtitle=paste0("",
       #caption=paste0("") 

# non_sp selection ----

non_sp <- rownames(sp.taxa.obj[sp.taxa.obj$SL != max(sp.taxa.obj$SL),])


if(identical(rownames(sp.taxa.obj), rownames(full.taxa.obj))) {
  sp_non_sp <- as_tibble(sp.taxa.obj[rownames(sp.taxa.obj) %in% non_sp, ])
  full_non_sp <- as_tibble(full.taxa.obj[rownames(full.taxa.obj) %in% non_sp, ])
}

#  quantitative comparison ----

n_na <- function(x) ( sum(is.na(x))) # sum(is.na(sp_non_sp[,2]))

na_number <- data.frame(
  sp = apply(sp_non_sp[ -9], 2, n_na),
  full = apply(full_non_sp[-9], 2, n_na))

na_number$Level <- rownames(na_number)
na_number$diff <- na_number$full - na_number$sp

na_number <- melt(na_number, id.vars = c('Level', 'diff'),
               variable.name = 'DataBase',
               value.name = 'n')

na_number$Level <- factor(na_number$Level, levels = TL)

ggplot(subset(na_number, diff != 0), aes(Level, n, fill = DataBase, group = DataBase)) + 
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  geom_text(data = subset(na_number, DataBase == 'sp' & diff != 0), aes(label=diff), size=4, vjust=1, color = 'black') +
  scale_fill_brewer(palette = "Set1") +
  labs(x="Rank", y="# ASVs (Undetermined | Unclassified)", 
       title = 'Number of NAs ASVs across full and species BOLD databases')

# count the level of change ----
x <- sp_non_sp$SL
y <- full_non_sp$SL

diff_x_y <- data.frame(
            Change = -6:4,# data.frame(table(x - y))[,1],
            Freq = data.frame(table(x - y))[,2],
            stringsAsFactors = FALSE)
# 

diff_x_y$Ref <- ''
diff_x_y[diff_x_y[,1] > 0, 'Ref'] <- 'full'
diff_x_y[diff_x_y[,1] < 0, 'Ref'] <- 'sp'

diff_x_y$Change <- abs(diff_x_y$Change)
diff_x_y$Change <- factor(diff_x_y$Change, levels = 0:6)


ggplot(subset(diff_x_y, Change != 0), aes(x=Change, y=log10(Freq), group=Ref, fill=Ref, color=Ref)) + 
  geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") +
  geom_text(data = subset(diff_x_y, Change != 0), aes(label=Freq), size=4, vjust=1, color = 'black') +
  scale_color_brewer(palette = "Set1") +
  labs(title="BOLD species - BOLD full db taxonomy Difference"
       #subtitle=paste0("",
       #caption=paste0("") 
  )


#if(full_non_sp[,1] == 'root') {full_non_sp <- full_non_sp[-1]}
#if(sp_non_sp[,1] == 'root') {sp_non_sp <- sp_non_sp[-1]}

# .... 1
x_ <- NULL
for (i in 1:nrow(sp.taxa.obj)) {
  rl <- sp.taxa.obj$SL[i] + 1
  x_[[i]] <- names(full.taxa.obj)[rl]
}

# 2.
for (i in 1:nrow(full.taxa.obj)) {
  rl <- full.taxa.obj$SL[i] + 1
  y_[[i]] <- names(full.taxa.obj)[rl]
}



# 1.
x_y_rank <- data.frame(ASV = rownames(sp.taxa.obj), 
                            sp= x_, full= y_, 
                       #SL_x = sp.taxa.obj$SL,
                       #SL_y = full.taxa.obj$SL,
                       x_y = sp.taxa.obj$SL - full.taxa.obj$SL,
                       stringsAsFactors = FALSE)


# 2.
x_y_rank_m <-  melt(x_y_rank, id.vars = c('ASV', 'x_y'),
                     variable.name = 'DataBase',
                     value.name = 'Rank')

x_y_rank_m$Rank <-  factor(x_y_rank_m$Rank, levels = TL)   
x_y_rank_m$DataBase <-  factor(x_y_rank_m$DataBase, levels = c('sp', 'full')) 
#aggregate(x_y_rank[,'DataBase'], by=list(x_y_rank[,'Rank']), FUN = table)

# Level change by rank
# agrupa los niveles por Rank y entonces ploteas el panel_grid

ggplot(subset(x_y_rank_m, x_y != 0), 
       aes(y=..count.., x=abs(x_y), fill = DataBase)) + 
  #geom_histogram(aes(y=..count..), position=position_dodge(), alpha=0.5, bins = 50) +
  geom_density(alpha=.7) + 
  scale_x_continuous(name = "Change", breaks  =  c(1:6)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x="change", y="Frequency of ASVs", 
       title = 'Level of Resolution between full (blue) and sp (red) BOLD dabases') +
  facet_wrap( ~ Rank , scales = 'free_y')

# alluvial

alluv <- subset(x_y_rank, x_y != 0 & full != 'root')
alluv$Ref <- 'sp'
alluv[alluv$x_y < 0, 'Ref'] <- 'full'
alluv$x_y <- abs(alluv$x_y)

alluv$sp <-  factor(alluv$sp, levels = TL[-1])
alluv$full <-  factor(alluv$full, levels = TL[-1]) 
alluv$Ref <- factor(alluv$Ref, levels = c('sp', 'full')) 

library(ggalluvial)

is_alluvia_form(alluv, axes = 1:3, silent = TRUE)

# library('ggalluvial')
ggplot(data = alluv,
       aes(axis1 = sp, axis2 = full, axis3 = Ref)) +
  geom_alluvium() + # aes(fill = Ref)
  geom_stratum() + 
  scale_fill_brewer(type = "qual", palette = "Set1") +
  geom_text(stat = "stratum", label.strata = TRUE, size = 3) +
  scale_x_discrete(limits = c("sp", "full", "Ref"), expand = c(.05, .05)) +
  # theme(legend.position = "bottom") +
  theme_minimal() +
  # use of lode controls
  geom_flow(aes(fill = sp, alpha = Ref), stat = "alluvium",
            aes.bind = TRUE, lode.guidance = "rightward")


# illustrate positioning
ggplot(data = alluv,
       aes(y = x_y, axis1 = sp, axis2 = full, axis3 = Ref, color = Ref)) +
  stat_stratum(geom = "errorbar") +
  geom_line(stat = "alluvium") +
  stat_alluvium(geom = "pointrange") +
  geom_text(stat = "stratum", label.strata = TRUE) +
  scale_x_discrete(limits = c("sp", "full", "Ref")) +
  theme_minimal()


