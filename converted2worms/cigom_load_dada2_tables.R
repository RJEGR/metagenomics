# cigom_load_dada2_tables.R
# Junio 2018
# Miguel Martinez | marmigues@gmail.com
#
# R functions that load parsed (cigom_parse_dada2_tables.R output) ASV/taxonomy tables from dada2
# load.mothurdp() - Carga el archivo 'wang'taxonomy' (sin bootstrap) generado en Mothur con la DB 'midori_unique_DB_0'
# load.mothurdp_boots() - Carga un archivo similar a load.mothurdp() pero que incluye los numeros de bootstrap. Ademas permite filtrar la asignacion por valor de bootstrap
# load.mothurdp_bworms() - Carga el archivo 'wang taxonomy' generado en Mothur con la DB 'ictio_coi_sanger114'. El archivo incluye el bootstrap y puede ser usado para filtrar la asignacion. 
# load.asvcounts() - Carga el archivo de conteos de ASVs con IDs codificados (ASV###)

# File for testing ---- 

# # Run and test
# runID <- #"run06_20170713_COI"
#     "run09_20180511_COI"
# run_n_test <- "run09_test08"
# # Taxonomy (assigned with mothur)
# mothurdp <- paste0("~/cigom/analisis_datos/ASV_analysis/", runID, "/", run_n_test, "/mothur_rdp_CLASSIFY/", run_n_test, "_asv.midori_unique_DB_0.wang.taxonomy")
# # coded-ASV counts
# asvcounts <- paste0("~/cigom/analisis_datos/ASV_analysis/", runID, "/", run_n_test, "/", run_n_test, "_asv_counts.tsv")
# Bootstrap cutoff
# boots <- 80

## Load mothur-rdp-classify ----

"load.mothurdp" <- function(mothurdp) {
    # Funcion que carga el archivo 'wang'taxonomy' generado en Mothur con la DB 'midori_unique_DB_0'
    # 'mothurdp' should be the complete path with the complete file name
    
    # Require
    require(tidyr)
    require(stringr)
    
    # Oficial taxonomy ranks
    TL <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Load
    rdp_table <- read.table(mothurdp, sep="\t", header=F, stringsAsFactors=F, na.strings=NA)
    names(rdp_table) <- c("ASV", "Taxonomy")
    
    ## Parse taxonomy file (separate each level into a column)
    # removes bootstrap values
    # removes last semicolon (of taxonomy)
    # 'root' is removed
    rdp_table$Taxonomy <- rdp_table$Taxonomy %>% str_replace_all(c("\\([0-9]*\\)"="", "\\([0-9]*\\.[0-9]*\\)"="", ";$"="", "^root;"=""))
    # Use tidyr to separate taxonomy
    rdp_df <- rdp_table %>% separate(Taxonomy, TL, sep=";", extra = "drop", fill = "right")
    
    # Check Kingdom is all Eukaryota
    #table(rdp_df$Kingdom)
    # Eukaryota 
    # 10843
    
    # Check how many are without assignation
    apply(rdp_df[,TL], 2, function(x) sum(is.na(x)))
    # Kingdom  Phylum   Class   Order  Family   Genus Species 
    # 0        0        0       0       0       0       0
    
    # Change 'unclassified' at Kingdom and Phylum level for "Undetermined"
    rdp_df$Phylum <- sub(".*nclassified.*", "Undetermined", rdp_df$Phylum, perl=TRUE)
    rdp_df$Kingdom[grep("Undetermined", rdp_df$Phylum, perl=TRUE)] <- "Undetermined"
        
    
    # Count how many undetermined (Phyla summary)
    #table(rdp_df$Phylum)
    # Annelida      Arthropoda    Chaetognatha        Chordata        Cnidaria      Ctenophora   Echinodermata        Mollusca Platyhelminthes    Undetermined 
    # 50            3012             469             301             380               1              73             404              16            6137    run09_test08
    # 13             671             125              48              59              NA              16              56               5            1279    run06_test02
    
    # Return rdp_df
    assign("rdp_df", rdp_df, .GlobalEnv)
    
} #ends load.mothurdp

## Load mothur-rdp-boots-classify ----

"load.mothurdp_boots" <- function(mothurdp, boots) {
    # Funcion que carga un archivo similar a load.mothurdp() pero que incluye los numeros de bootstrap. Ademas permite filtrar la asignacion por valor de bootstrap
    # 'mothurdp' should be the complete path with the complete file name
    # 'boots' is the bootstrap cutoff
    
    # Require
    require(tidyr)
    require(stringr)
    library(tidyverse)
    
    # Oficial taxonomy ranks
    TL <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Load
    rdp_table <- read.table(mothurdp, sep="\t", header=F, stringsAsFactors=F, na.strings=NA)
    names(rdp_table) <- c("ASV", "Taxonomy")
    
    ## Parse taxonomy file (extract bootstrap values)
    # separate taxonomy levels
    boots_table <- rdp_table %>% separate(Taxonomy, paste0(c("root", TL), "_boots"), sep=";", extra = "drop", fill = "right")
    # remove 'root'
    boots_table$root_boots <- NULL
    # keep only numbers
    boots_table <- boots_table %>% mutate_all(funs(str_replace_all(., c(".+\\(" = "" , "\\)$"=""))))
    # save IDs for later and remove
    boots_rows <- boots_table$ASV
    boots_table$ASV <- NULL
    # Ensure we have numerical columns (convert to matrix)
    boots_table <- apply(boots_table, 2, as.numeric)
    # FILTER (step 1/2) values below bootstrap threshold (convert them to NA)
    boots_table[boots_table < boots] <- NA
    # Format filtered DF
    boots_table <- data.frame(boots_table)
    rownames(boots_table) <- boots_rows
    head(boots_table)
    
    ## Parse taxonomy file (separate each level into a column)
    # removes bootstrap values
    # removes last semicolon (of taxonomy)
    # 'root' is removed
    rdp_table$Taxonomy <- rdp_table$Taxonomy %>% str_replace_all(c("\\([0-9]*\\)"="", "\\([0-9]*\\.[0-9]*\\)"="", ";$"="", "^root;"=""))
    # Use tidyr to separate taxonomy
    rdp_df <- rdp_table %>% separate(Taxonomy, TL, sep=";", extra = "drop", fill = "right")
    head(rdp_df)
    
    ## Combine DFs
    mypaste <- function(x,y) paste(x, y, sep="_")
    rdp_boots <- mapply(mypaste, rdp_df[,-1], boots_table)
    
    ## FILTER (step 2/2) by bootstrap, use NA's appended from numeric columns
    rdp_boots <- data.frame(rdp_boots) %>% mutate_all(funs(str_replace_all(., c(".+_NA"= NA, "_\\d+.+$"=""))))
    
    ## Organize DF
    rdp_boots$ASV <- rdp_df$ASV
    rdp_boots <- rdp_boots[,c("ASV", TL)]
    
    # Percentaje of ASVs with assignation (per level)
    apply(rdp_boots[,TL], 2, function(x) round((nrow(rdp_table) - sum(is.na(x)))/nrow(rdp_table)*100, 2))
    # Kingdom  Phylum   Class   Order  Family   Genus Species 
    # 100.00    42.21   37.55   36.84   35.21   33.28  ## boots > 99
    # 100.00    63.17   48.13   44.57   43.35   41.18   39.75  ## boots > 80
    
    ## Change unclassfied to undetermined
    rdp_boots$Phylum <- sub(".*nclassified.*", "Undetermined", rdp_boots$Phylum, perl=TRUE)
    rdp_boots$Kingdom[grep("Undetermined", rdp_boots$Phylum, perl=TRUE)] <- "Undetermined"
    
    # Count how many undetermined (Phyla summary)
    #nrow(rdp_boots[which(rdp_boots$Phylum %in% "Undetermined"), ])
    ## 0     #run012

    # Clean space
    rm(mothurdp)
    
    # Return rdp_boots
    assign("rdp_df", rdp_boots, .GlobalEnv)
    
} #ends load.mothurdp_boots


## Load mothur-rdp-boots-worms-classify ----

"load.mothurdp_bworms" <- function(mothurdp, boots) {
    # Funcion que carga el archivo 'wang taxonomy' generado en Mothur con la DB 'ictio_coi_sanger114'. El archivo incluye el bootstrap y puede ser usado para filtrar la asignacion. 
    # 'mothurdp' should be the complete path with the complete file name
    # 'boots' is the bootstrap cutoff
    
    # Require
    require(tidyr)
    require(stringr)
    library(tidyverse)
    
    # Oficial taxonomy ranks
    TL <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Load
    rdp_table <- read.table(mothurdp, sep="\t", header=F, stringsAsFactors=F, na.strings=NA)
    names(rdp_table) <- c("ASV", "Taxonomy")
    
    ## Parse taxonomy file (extract bootstrap values)
    # separate taxonomy levels
    boots_table <- rdp_table %>% separate(Taxonomy, paste0(TL, "_boots"), sep=";", extra = "drop", fill = "right")
    # keep only numbers
    boots_table <- boots_table %>% mutate_all(funs(str_replace_all(., c(".+\\(" = "" , "\\)$"=""))))
    # save IDs for later and remove
    boots_rows <- boots_table$ASV
    boots_table$ASV <- NULL
    # Ensure we have numerical columns (convert to matrix)
    boots_table <- apply(boots_table, 2, as.numeric)
    # FILTER (step 1/2) values below bootstrap threshold (convert them to NA)
    boots_table[boots_table < boots] <- NA
    # Format filtered DF
    boots_table <- data.frame(boots_table)
    rownames(boots_table) <- boots_rows
    head(boots_table)
    
    ## Parse taxonomy file (separate each level into a column)
    # removes bootstrap values
    # removes last semicolon (of taxonomy)
    rdp_table$Taxonomy <- rdp_table$Taxonomy %>% str_replace_all(c("\\([0-9]*\\)"="", "\\([0-9]*\\.[0-9]*\\)"="", ";$"=""))
    # Use tidyr to separate taxonomy
    rdp_df <- rdp_table %>% separate(Taxonomy, TL, sep=";", extra = "drop", fill = "right")
    head(rdp_df)
    
    ## Combine DFs
    mypaste <- function(x,y) paste(x, y, sep="_")
    rdp_boots <- mapply(mypaste, rdp_df[,-1], boots_table)
    
    ## FILTER (step 2/2) by bootstrap, use NA's appended from numeric columns
    rdp_boots <- data.frame(rdp_boots) %>% mutate_all(funs(str_replace_all(., c(".+_NA"= NA, "_\\d+.+$"=""))))
    
    ## Organize DF
    rdp_boots$ASV <- rdp_df$ASV
    rdp_boots <- rdp_boots[,c("ASV", TL)]
    
    # Percentaje of ASVs with assignation (per level)
    # apply(rdp_boots[,TL], 2, function(x) round((nrow(rdp_table) - sum(is.na(x)))/nrow(rdp_table)*100, 2))

    ## Change unclassfied to undetermined
    rdp_boots$Phylum <- sub(".*nclassified.*", "Undetermined", rdp_boots$Phylum, perl=TRUE)
    rdp_boots$Kingdom[grep("Undetermined", rdp_boots$Phylum, perl=TRUE)] <- "Undetermined"
    
    # Count how many undetermined (Phyla summary)
    #nrow(rdp_boots[which(rdp_boots$Phylum %in% "Undetermined"), ])
    ## 0     #run012
    
    # Clean space
    rm(mothurdp)
    
    # Return rdp_boots
    assign("rdp_df", rdp_boots, .GlobalEnv)
    
} #ends load.mothurdp_bworms


## Load coded-ASVs ----

"load.asvcounts" <- function(asvcounts) {
    # Funcion que carga el archivo de conteos de ASVs con IDs codificados (ASV###)
    # 'asvcounts' should be the complete path with the complete file name
    
    # Require
    require(readr)
    
    # Load with readr (faster)
    asv_table <- read_tsv(asvcounts, col_names = TRUE)
    # Load with base-R (slow)
    #asv_table <- read.table(asvcounts, sep="\t", header=T, stringsAsFactors=F, na.strings=NA)
    dim(asv_table)
    # 53 10844  run09_test08
    
    # Return table
    assign("asv_table", asv_table, .GlobalEnv)

} #Ends load.asvcounts

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
} #Ends "read.wangtax"