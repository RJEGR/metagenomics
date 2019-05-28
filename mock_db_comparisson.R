#!/usr/bin/env Rscript

rm(list=ls()); 

# Functions ----
# Back the SL position per database
sl_back <- function(x) {
  x_ <- NULL
  for (i in 1:nrow(x)) {
  rl <- x$SL[i] + 1
  x_[[i]] <- list(rank=names(x)[rl], linage=x[i,rl]) }
  x_ <- do.call(rbind, x_)
  return(x_)
}

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


path <- '/Users/cigom/metagenomics/sanger_assign_db'

setwd(path)

bold_fullf <- list.files(pattern = "*ALL.wang.taxonomy")
bold_spf <- list.files(pattern = "*BOLD_public_species.wang.taxonomy")
gbank_f <- list.files(pattern = "*Coarbitrator_COI_nuc_curated.wang.taxonomy")
midori_f <- list.files(pattern = "*midori_unique_DB_0.wang.taxonomy")

# Load files ----
bold_full <- read_rdp(bold_fullf)
bold_sp <- read_rdp(bold_spf)
gbank <- read_rdp(gbank_f)
midori <- read_rdp(midori_f)

TL <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#colnames(full.taxa.obj) <- c(TL, 'SL')

# Level Resolution vis ----

SL <- rbind(data.frame(SL = bold_full$SL,
                       DataBase = 'full'),
            data.frame(SL = bold_sp$SL,
                       DataBase = 'sp'),
            data.frame(SL = gbank$SL,
                       DataBase = 'gbank'),
            data.frame(SL = midori$SL,
                       DataBase = 'midori')
            )


# SL <- SL[SL$SL > 0 ,]
SL$SL <- factor(SL$SL, levels = names(table(SL$SL)))

# and
SL_agg <- aggregate(SL[,'DataBase'], by=list(SL[,'SL']), FUN = table)
SL_agg <- data.frame(level = SL_agg[,1], SL_agg[,2])
SL_agg$Rank <- c(TL[-1], rep('Extra', 6))
SL_aggM <- melt(SL_agg)
SL_aggM$Rank <- factor(SL_aggM$Rank, levels = c(TL[-1], 'Extra'))

levels <- c("full","sp", "gbank", "midori")

SL_aggM$variable <- factor(SL_aggM$variable, levels = levels)

ggplot(SL_aggM, aes(x=Rank, y=value, group=variable, fill=variable, color=variable)) + 
  geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") +
  geom_text(data = SL_aggM, aes(label=value), size=4, vjust=1, color = 'black') +
  scale_color_brewer(palette = "Set1") +
  labs(title="db taxonomy Difference",
       x = 'Taxonomy level',
       y = '# ASVs')





#subtitle=paste0("",
#caption=paste0("") 

# Prepare alluvial ----
# .... 1
x_ <- NULL
for (i in 1:nrow(sp.taxa.obj)) {
  rl <- sp.taxa.obj$SL[i] + 1
  x_[[i]] <- names(full.taxa.obj)[rl]
}

# 2.
y_ <- NULL
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
x_y_rank_m$DataBase <-  factor(x_y_rank_m$DataBase, levels = levels)

#aggregate(x_y_rank[,'DataBase'], by=list(x_y_rank[,'Rank']), FUN = table)

# Level change by rank ----
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

# alluvial ----

alluv <- subset(x_y_rank, x_y != 0 & full != 'root')
alluv$Ref <- 'full'
alluv[alluv$x_y < 0, 'Ref'] <- 'sp'
alluv$x_y <- abs(alluv$x_y)

alluv$sp <-  factor(alluv$sp, levels = TL[-1])
alluv$full <-  factor(alluv$full, levels = TL[-1]) 
alluv$Ref <- factor(alluv$Ref, levels = levels) 

library(ggalluvial)

is_alluvia_form(alluv, axes = 1:3, silent = TRUE)

alluv <- alluv[order(alluv$sp, decreasing = TRUE), ]


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



#  quantitative comparison ----

n_na <- function(x) ( sum(is.na(x))) # sum(is.na(sp_non_sp[,2]))

bold_full <- read_rdp(bold_fullf)
bold_sp <- read_rdp(bold_spf) # se esta mandando esta version con el nivel root;Eukaryota
gbank <- read_rdp(gbank_f)
midori <- read_rdp(midori_f)

na_number <- list(
  sp = apply(bold_full[ -ncol(bold_full)], 2, n_na),
  full = apply(bold_sp[-ncol(bold_sp)], 2, n_na),
  gbank = apply(gbank[-ncol(gbank)], 2, n_na),
  midori = apply(midori[-ncol(midori)], 2, n_na))

na_number$Level <- rownames(na_number)
na_number$diff <- na_number$full - na_number$sp
na_number_m <- melt(na_number, id.vars = c('Level', 'diff'),
                    variable.name = 'DataBase',
                    value.name = 'n')

na_number_m$Level <- factor(na_number_m$Level, levels = TL)
na_number_m$DataBase <- factor(na_number_m$DataBase, levels = levels)

ggplot(subset(na_number_m, diff != 0), aes(Level, n, fill = DataBase, group = DataBase)) + 
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  geom_text(data = subset(na_number_m, DataBase == 'sp' & diff != 0), aes(label=abs(diff)), size=4, vjust=1, color = 'black') +
  scale_fill_brewer(palette = "Set1") +
  labs(x="Rank", y="# ASVs (Undetermined | Unclassified)", 
       title = 'Number of NAs ASVs across full and species BOLD databases')

# count the level of change ----
x <- sp.taxa.obj$SL
y <- full.taxa.obj$SL

diff_x_y <- data.frame(
  Change = -7:6,# data.frame(table(x - y))[,1],
  Freq = data.frame(table(x - y))[,2],
  stringsAsFactors = FALSE)
# 

diff_x_y$Ref <- ''
diff_x_y[diff_x_y[,1] > 0, 'Ref'] <- 'sp'
diff_x_y[diff_x_y[,1] < 0, 'Ref'] <- 'full'

diff_x_y$Change <- abs(diff_x_y$Change)
diff_x_y$Change <- factor(diff_x_y$Change, levels = 0:7)

diff_x_y$Ref <- factor(diff_x_y$Ref, levels = levels)

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

