read_rdp <- function(taxonomy.file) {
  .cran_packages <- c('stringr', 'dplyr')
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
  
  # if(tax[,1] == 'root') {
  #   tax <- tax[-1]
  #   # ROW/COL names definition
  #   rank.names <- vector(max.rank, mode="character")
  #   for (i in 1:max.rank) { rank.names[i] <- paste("Rank", i, sep="_") }
  #   colnames(tax) <- rank.names
  #   } else {
  # ROW/COL names definition

  rank.names <- vector(max.rank, mode="character")
  for (i in 1:max.rank) { rank.names[i] <- paste("Rank", i, sep="_") }
  colnames(tax) <- rank.names
  
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
  
  # tax$SL_rank <- 'rank'
  # for (i in 1:nrow(tax)) {
  #   rl <- tax$SL[i] + 1
  #   tax$SL_rank[i] <- names(tax)[rl] }
  
  return(tax) 
  
  
}