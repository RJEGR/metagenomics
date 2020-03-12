read_rdp <- function(taxonomy.file, header = FALSE) {
  .cran_packages <- c('stringr', 'dplyr')
  sapply(c(.cran_packages), require, character.only = TRUE, quietly = TRUE)
  taxonomy.obj <- read.csv(taxonomy.file, header=header, sep="\t", stringsAsFactors=FALSE)
  tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
  
  # Using the max rank assignation to names the taxonomy object
  max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
  taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
  taxonomy <- as.data.frame(t(taxonomy))
  tax <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)
  
  # Remove unclassified and replace underscores by spaces
  tax <- mutate_all(data.frame(tax), funs(str_replace_all(., c(".+_unclassified"=NA_character_, "_"=" ", "Unclassified"=NA_character_, "NA"=NA_character_, "sp."=NA_character_, "Undetermined"=NA_character_, "Unknown"=NA_character_))))
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
  # Tag starting level (Resolution) (aka calculate resolution)
  tax$Resolution <- max.rank
  
  for(t in max.rank:2){
    rr <- t-1  #rank real 
    rc <- t+1  #column-index (rank) para corregir
    if(t == max.rank){
      tax$Resolution[is.na(tax[,rc])] <- rr
    } else {
      tax$Resolution[is.na(tax[,rc]) & tax$Resolution <= t ] <- rr
    }
  }
  
  # tax$Resolution_rank <- 'rank'
  # for (i in 1:nrow(tax)) {
  #   rl <- tax$Resolution[i] + 1
  #   tax$Resolution_rank[i] <- names(tax)[rl] }
  
  return(tax) 
  
  
}


# get boots ----
# cat("\n 2. Bootstrap stats ... \n")
boots_rdp <- function(taxonomy.file) {
  
  boots.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
  tax.split <- strsplit(boots.obj[, ncol(boots.obj)], ";")
  max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
  
  boots0 <- sapply(tax.split, "[", c(1:max.rank)) # Using the max rank assignation to names the taxonomy object
  boots0 <- as.data.frame(t(boots0))
  
  boots <- as.data.frame(apply(boots0, 2, function(x) gsub("[A-z||()]", "",  x, perl=TRUE)), stringsAsFactors = F)
  #boots <- apply(boots, 2, as.numeric)
  #boots <- data.frame(boots)
  #boots.obj[which(is.na(boots$V9)), 2]
  #boots[is.na(boots)] <- 0
  
  rownames(boots) <- boots.obj[,1]
  
  return(boots)
  
  
  # boots.backup <- boots
  
  # # Recorte de algun rank con asignacion root;
  # boots <- boots[apply(boots, 2, mean) != 100]
  # boots <- boots[apply(boots, 1, max) <= 100,] #remove rows redundant bootstrap
  # boots <- boots[apply(boots, 1, min) > 0,] # ""
  
}


# bold-sp vs bold-full ----
# bbold(filter(x, full == TL[2:6]), fasta_file = fasta_file, count_tbl = count_tbl)
bbold <- function(y, fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE) {
  
  require(Biostrings)
  
  x <- y[1]
  db_subset <- x[!duplicated(x)]
  # 1) abundance ----
  
  count.tbl0 <- read.table(count_tbl, row.names = 1)
  total <- sum(rowSums(count.tbl0))
  nreads <- rowSums(count.tbl0)
  
  # count.tbl <- data.frame(ASV = rownames(count.tbl0[rownames(count.tbl0) %in% db_subset, ]),
  #                         abund = rowSums(count.tbl0[rownames(count.tbl0) %in% db_subset, ]))
  
  if(rel_ab) {
    rabund <- data.frame(ASV = names(nreads), rabund = (nreads / total) * 100)
    
    rabund0 <- rabund[rabund$ASV %in% db_subset, ]
    rabund0 <- rabund0[match(db_subset, rownames(rabund0)),] 
    } else {
    rabund0 <- data.frame(ASV = names(nreads), rabund = nreads)
    rabund0 <- rabund0[match(db_subset, rownames(rabund0)),] 
    
  }
  
  if(identical(rownames(rabund0), db_subset)) {
    
    count.tbl <- data.frame(ASV = rownames(count.tbl0)[which(rownames(count.tbl0) %in% db_subset)],
                            abund = rabund0$rabund)
    count.tbl <- count.tbl[match(db_subset, count.tbl$ASV),]
    }
  
# sanity check
  if (identical(as.character(count.tbl$ASV), db_subset) != TRUE) stop('... The order in the count.tbl is distinct than present subset and is needed... Stopping. 1');
  
  # 2) sequence size ----
  seqs0 <- readDNAStringSet(fasta_file)
  seqs <- seqs0[names(seqs0) %in% db_subset]
  seqs <- seqs[match(db_subset, names(seqs)),]

# sanity check
  if (identical(names(seqs), db_subset) == FALSE) stop('... Stopping. 2');
  
  # 3) taxonomy db distribution ----

# parse data.frame ----

  seqs <- seqs[match(db_subset, names(seqs)),]
  #count.tbl <- count.tbl[match(db_subset, count.tbl$ASV),]


  if(identical(names(seqs), as.character(count.tbl$ASV))) {
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
  names(full_set) <- paste("incomplete", TL2, sep="_")
  names(sp_set) <- paste("complete", TL2, sep="_")


# sanity check


  if(identical(rownames(full_set), rownames(sp_set))) {
  if(identical(rownames(full_set), db_subset)) {
    if(identical(rownames(full_set), out$ASV)) {
      out <- cbind(out, full_set, sp_set)
      out <- out[order(out$abund, decreasing = FALSE),]
      }
     }
    }
  
  # return(out)
  
  # additional ----
  
  
  if(exists("x_y_rank")) {
    require(dplyr)
    out %>%
      group_by(ASV) %>%
      inner_join(x_y_rank, by = 'ASV') %>%
      select(-complete, -incomplete) %>%
      db_color() -> out0
    return(out0)
  } else {
    out0 <- out
    return(out0)
  }
  
}

# alluvial colored ----
db_color <- function(x) {
  # use for object alluv <- subset(x_y_rank, x_y != 0 & full != 'root')
  x$Ref <- 'incomplete'
  x[x$x_y > 0, 'Ref'] <- 'complete' 
  return(x)
}

# count out levels ----
aglom_ab <- function(x,rank) {
  tax_g <- aggregate(x[,'abund'],  by = list(x[, rank]), FUN = sum)
  n_asvs <- aggregate(x[,rank],  by = list(x[, rank]), FUN = length)
  min_seq_size <- aggregate(x[,'seq_size'],  by = list(x[, rank]), FUN = min)
  max_nsamples <- aggregate(x[,'nsamples'],  by = list(x[, rank]), FUN = max)
  
  if(identical(tax_g[,1], n_asvs[,1])){
    if(identical(tax_g[,1], min_seq_size[,1])){
      tax_out <- data.frame(lineage = tax_g[,1], Abundance = tax_g[,2], Size = n_asvs[,2] , mseq_size = min_seq_size[,2], Mnsamples = max_nsamples[,2])
    }
  }
  # tax_out <- data.frame(lineage = tax_g[,1], Size = tax_g[,2])
  return(tax_out)
}

# agglomerate nsamples per rank and add children level to rank (REMOVE)

# get back the last rank based on the Resolution ----

lrank <- function(x) {
  x_ <- NULL
  for (i in 1:nrow(x)) {
    rl <- x$Resolution[i] + 1
    x_[[i]] <- names(x)[rl]
  }
  return(x_)
}

# get back the last lineage based on the Resolution  ----

llineage <- function(x) {
  x_ <- NULL
  for (i in 1:nrow(x)) {
    rl <- x$Resolution[i] + 1
    # x_[[i]] <- list(rank=names(x)[rl], linage=x[i,rl]) }
    x_[[i]] <- x[i,rl]
  }
  return(x_)
}

# additional step to bbold ---
bbold_old <- function(x, x_y_rank) {
  x %>%
    group_by(ASV) %>%
    inner_join(x_y_rank, by = 'ASV') %>%
    select(-complete, -incomplete) %>%
    db_color() -> R_out
  return(R_out)
}

# Rank levels to numeric


rank_n <- function(x) {
  x %>%
    mutate(Rank_n = if_else(Rank == 'Domain', "1",
                            if_else(Rank == 'Kingdom', "2",
                                    if_else(Rank == 'Phylum', "3",
                                            if_else(Rank == 'Class', "4",
                                                    if_else(Rank == 'Order', "5",
                                                            if_else(Rank == 'Family', "6",
                                                                    if_else(Rank == 'Genus', "7",
                                                                            if_else(Rank == 'Species', "8", "")
                                                                    )))))))) -> out
  out$Rank_n <- as.numeric(out$Rank_n)
  out <- out %>% group_by(ASV)
  return(out)
}


# bbold single version

bbold_ <- function(y, fasta_file = fasta_file, count_tbl = count_tbl, rel_ab = TRUE) {
  
  options(stringsAsFactors = FALSE)
  require(Biostrings)
  
  x <- y[,1]
  db_subset <- as.character(x[!duplicated(x)])
  # 1) abundance ----
  
  count.tbl0 <- read.table(count_tbl)
  total <- sum(rowSums(count.tbl0))
  nreads <- rowSums(count.tbl0)
  
  nsamples <- data.frame(nsamples = c(apply(count.tbl0, 1, function(x) {sum(x > 0)})))

  if(rel_ab) {
    rabund <- data.frame(ASV = names(nreads), rabund = (nreads / total) * 100)
    
    rabund0 <- rabund[rabund$ASV %in% db_subset, ]
    rabund0 <- rabund0[match(db_subset, rownames(rabund0)),] 
  } else {
    rabund0 <- data.frame(ASV = names(nreads), rabund = nreads)
    rabund0 <- rabund0[match(db_subset, rownames(rabund0)),] 
    
  }
  
  if(identical(rownames(rabund0), db_subset)) {
    if(identical(rownames(nsamples), db_subset)) {
      
      count.tbl <- data.frame(ASV = rownames(count.tbl0)[which(rownames(count.tbl0) %in% db_subset)], abund = rabund0$rabund,  nsamples)
      
      count.tbl <- count.tbl[match(db_subset, count.tbl$ASV),]
      
    }
  }
  
  # sanity check
  if (identical(as.character(count.tbl$ASV), db_subset) != TRUE) stop('... The order in the count.tbl is distinct than present subset and is needed... Stopping. 1');
  
  # 2) sequence size ----
  seqs0 <- readDNAStringSet(fasta_file)
  seqs <- seqs0[names(seqs0) %in% db_subset]
  seqs <- seqs[match(db_subset, names(seqs)),]
  
  # sanity check
  if (identical(names(seqs), db_subset) == FALSE) stop('... Stopping. 2');
  
  # 3) taxonomy db distribution ----
  
  # parse data.frame ----
  
  seqs <- seqs[match(db_subset, names(seqs)),]
  #count.tbl <- count.tbl[match(db_subset, count.tbl$ASV),]
  
  
  if(identical(names(seqs), as.character(count.tbl$ASV))) {
    out <- data.frame(ASV = names(seqs),
                      seq_size = width(seqs),
                      abund = count.tbl$abund,
                      nsamples = count.tbl$nsamples,
                      select(tax, -ASV),
                      stringsAsFactors = FALSE)
  }
  
   return(out)
  
}
