read_rdp <- function(taxonomy.file) {
  .cran_packages <- c('stringr', 'dplyr')
  sapply(c(.cran_packages), require, character.only = TRUE, quietly = TRUE)
  taxonomy.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
  
  # Using the max rank assignation to names the taxonomy object
  max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
  taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
  taxonomy <- as.data.frame(t(taxonomy))
  tax <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)
  
  # Remove unclassified and replace underscores by spaces
  tax <- mutate_all(data.frame(tax), funs(str_replace_all(., c(".+_unclassified"=NA_character_, "_"=" ", "Unclassified"=NA_character_, "NA"=NA_character_))))
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

# bold-sp vs bold-full ----
# bbold(filter(x, full == TL[2:6]), fasta_file = fasta_file, count_tbl = count_tbl)
bbold <- function(y, fasta_file = fasta_file, count_tbl = count_tbl) {
  
  require(Biostrings)
  
  x <- y$ASV
  db_subset <- x[!duplicated(x)]
  # 1) abundance ----
  
  count.tbl0 <- read.table(paste0(path_BOLD, "/",count_tbl))

  # count.tbl <- data.frame(ASV = rownames(count.tbl0[rownames(count.tbl0) %in% db_subset, ]),
  #                         abund = rowSums(count.tbl0[rownames(count.tbl0) %in% db_subset, ]))
  
  total <- sum(rowSums(count.tbl0))
  nreads <- rowSums(count.tbl0)
  # dim(out <- out[out$abund > 1, ]) 
  rabund <- data.frame(ASV = names(nreads),
                       rabund = (nreads / total) * 100)
  
  
  rabund0 <- rabund[rabund$ASV %in% db_subset, ]
  rabund0 <- rabund0[match(db_subset, rownames(rabund0)),]
  
  if(identical(rownames(rabund0), db_subset)) {
    
    count.tbl <- data.frame(ASV = rownames(count.tbl0[rownames(count.tbl0) %in% db_subset, ]),
                            abund = rabund0$rabund) }
  
  count.tbl <- count.tbl[match(db_subset, count.tbl$ASV),]
  
# sanity check
  if (identical(count.tbl$ASV, db_subset) == FALSE) stop('... Stopping. 1');
  
  # 2) sequence size ----
  seqs0 <- readDNAStringSet(paste0(path_BOLD,'/',fasta_file))
  seqs <- seqs0[names(seqs0) %in% db_subset]
  seqs <- seqs[match(db_subset, names(seqs)),]

# sanity check
  if (identical(names(seqs), db_subset) == FALSE) stop('... Stopping. 2');
  
  # 3) taxonomy db distribution ----

# parse data.frame ----

  seqs <- seqs[match(db_subset, names(seqs)),]
  #count.tbl <- count.tbl[match(db_subset, count.tbl$ASV),]


  if(identical(names(seqs), count.tbl$ASV)) {
  out <- data.frame(ASV = names(seqs),
                    seq_size = width(seqs),
                    abund = count.tbl$abund,
                    stringsAsFactors = FALSE)
  }

# summary(out$abund)


  # get the lineage assignation per db ----


  full_set <- select(full.taxa.obj[rownames(full.taxa.obj) %in% db_subset,], TL[-1])
  sp_set <- select(sp.taxa.obj[rownames(sp.taxa.obj) %in% db_subset,], TL[-1])


  full_set <- full_set[match(db_subset, rownames(full_set)),]

  sp_set <- sp_set[match(db_subset, rownames(sp_set)),]


  TL2 <- c("D", "K", "P", "C", "O", "F", "G", "S")
  names(full_set) <- paste("full", TL2, sep="_")
  names(sp_set) <- paste("sp", TL2, sep="_")


# sanity check


  if(identical(rownames(full_set), rownames(sp_set))) {
  if(identical(rownames(full_set), db_subset)) {
    if(identical(rownames(full_set), out$ASV)) {
      out <- cbind(out, full_set, sp_set)
      out <- out[order(out$abund, decreasing = FALSE),]
      }
     }
    }
  
  return(out)
  
  
}
