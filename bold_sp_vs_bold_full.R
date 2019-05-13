#!/usr/bin/env Rscript

rm(list=ls()); 

# Functions ----

read_rdp <- function(taxonomy.file) {
  .cran_packages <- c('dplyr', 'stringr')
  sapply(c(.cran_packages), require, character.only = TRUE)
  taxonomy.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
 
  # Using the max rank assignation to names the taxonomy object
  max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
  taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
  taxonomy <- as.data.frame(t(taxonomy))
  tax <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)
  
  # Remove unclassified and replace underscores by spaces
  tax <- mutate_all(data.frame(tax), funs(str_replace_all(., c(".+_unclassified"=NA, "_"=" "))))
  rownames(tax) <- taxonomy.obj[,1]
  # colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  max.rank <-  ncol(tax) -1 
  # Tag starting level (SL) (aka calculate resolution)
  tax$SL <- max.rank
  
  for(t in max.rank:2){
    rr <- t-1  #rank real 
    rc <- t+1  #column-index (rank) para corregir
    if(t == max.rank){
      tax$SL[is.na(tax[,rc])] <- rr
    } else {
      tax$SL[is.na(tax[,rc]) & tax$SL <= t ] <- rr
    }
  }
  
  return(tax) 
  

}


# Set filenames ----
path_BOLD <- '/Users/cigom/metagenomics/COI/species_resolution_per_db'
bold_all <- 'run014_t2_ASVs.ALL.wang.taxonomy'
bold_sp <- 'run014_t2_ASVs.BOLD_public_species_v02.wang.taxonomy'


# Load files ----
full.taxa.obj <- read_rdp(paste0(path_BOLD,'/',bold_all))
sp.taxa.obj <- read_rdp(paste0(path_BOLD,'/',bold_sp))

# Special procedure to full db ----
full.taxa.obj$SL <- full.taxa.obj$SL - 1
full.taxa.obj <- full.taxa.obj[,c(1,3:10)]

# Level Resolution vis ----

SL <- rbind(data.frame(SL = full.taxa.obj$SL,
                 DataBase = 'full'),
            data.frame(SL = sp.taxa.obj$SL,
                       DataBase = 'SP'))

ggplot(SL, aes(y=..count.., x=SL, fill = DataBase)) + 
  geom_density(alpha=.7) + scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(x="Level of Resolution", y="Frequency of ASVs", 
       title = 'Level of Resolution between full (blue) and sp (red) BOLD dabases')

# 
TL <- c("root","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

colnames(full.taxa.obj) <- c(TL, 'SL')
colnames(sp.taxa.obj) <- c(TL, 'SL')

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


# done