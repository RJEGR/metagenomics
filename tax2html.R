
# 
# Junio 2019
# Ricardo Gomez
# Rscript --vanilla ~/Documents/GitHub/metagenomics/tax2html.R wang.taxonomy count_tbl fasta_file

rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors = FALSE)

wang.taxonomy <- args[1]
count_tbl <- args[2]
fasta_file <- args[3]

wang.taxonomy <- 'run014_t2_ASVs.ALL.wang.taxonomy'
count_tbl <- 'run014_t2_ASVs_count.table'
fasta_file <- 'run014_t2_ASVs.fasta'

# ------ functions

# Include the convertion to retrive BOLD id or BOLD BIND with

boldid_ <- function(x) {
  require(taxize)
  boldid <- get_boldid(x, verbose = FALSE, check = FALSE)
  id <- paste0('<a href=', attributes(boldid)$uri, '>', boldid[1], '</a>')
  return(id)
}

# 

boldid_2 <- function(x) {
  x_ <- NULL
  for (i in 1:nrow(x)) {
    rl <- x$SL[i] + 1
    x_[[i]] <- boldid_(x[i, rl])
  }
  return(x_)
}

# Path ----

# path <- getwd()

path <- '/Users/cigom/metagenomics/COI/species_resolution_per_db'

# Load files ----

tax.file <- list.files(path = path, full.names = TRUE, pattern = wang.taxonomy)
count.file <- list.files(path = path, full.names = TRUE, pattern = count_tbl)
fasta.file <- list.files(path = path, full.names = TRUE, pattern = fasta_file)


cat('File to process: ', tax.file, "\n")
cat('with ', count_tbl, 'and ', fasta_file)

TL <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# taxonomy.obj <- read.csv(tax.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
# tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
# 
# # Using the max rank assignation to names the taxonomy object
# 
# max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
# taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
# taxonomy <- as.data.frame(t(taxonomy))
# tax <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)
# 
# # Remove unclassified and replace underscores by spaces
# tax <- mutate_all(data.frame(tax), funs(str_replace_all(., c(".+_unclassified"=NA_character_, "_"=" ", "Unclassified"=NA_character_, "NA"=NA_character_))))
# rownames(tax) <- taxonomy.obj[,1]

# rank.names <- vector(max.rank, mode="character")
# for (i in 1:max.rank) { rank.names[i] <- paste("Rank", i, sep="_") }
# colnames(tax) <- rank.names

# set the last resolved lineage based on the SL

# Tag starting level (SL) (aka calculate resolution)
# max.rank <-  ncol(tax) -1  # ignorate 'root' label
# tax$SL <- max.rank
# 
# for(t in max.rank:2){
#   rr <- t-1  #rank real 
#   rc <- t+1  #column-index (rank) para corregir
#   if(t == max.rank){
#     tax$SL[is.na(tax[,rc])] <- rr
#   } else {
#     tax$SL[is.na(tax[,rc]) & tax$SL <= t ] <- rr
#   }
# }

# try
# boldid_(tax[1, max.rank])

# ----- 

source(file = "~/Documents/GitHub/metagenomics/readtx.R")

tax <- read_rdp(tax.file)

colnames(tax) <- c(TL, 'SL')

cat('Sample size is:', nrow(tax), "\n")

bold_tax <- boldid_2(tax)

out_tbl <- data.frame(select(tax, TL[-1]), boldid = bold_tax)

require(DT)

widget <- datatable(
  out_tbl, 
  escape=1,
  extensions = 'Buttons', options = list(
    pageLength = 25,  
    dom = 'Bfrtip',
    buttons = 
      list('copy', 'print', list(
        extend = 'collection',
        buttons = c('csv', 'excel'),
        text = 'Download'
      ))
    
  )
)

htmlwidgets::saveWidget(widget, paste0(tax.file, ".html"))

cat("\nDataTable conversion was done:", wang.taxonomy,"\n")

quit(save = "no")

# YOU CAN ALSO SUMMARY LINEAGE ITH NUMBER OF ASVS WITH THIS DISPOSAL

tax <- data.frame(ASV = rownames(tax), tax)

out0 <- bbold_(tax, fasta.file, count.file,  rel_ab = FALSE)


str(rank_out <- data.frame(aglom_ab(out0, 'Species')))


boldid_3 <- function(x) {
  x_ <- NULL
  for (i in 1:nrow(x)) {
    x_[[i]] <- boldid_(x$lineage[i])
  }
  return(x_)
}

# connect to a secured network (or error with Web authent)

boldid_(rank_out$lineage[2])


bold3_tax <- boldid_3(rank_out) # time demand!
out_tbl <- data.frame(rank_out, boldid = bold3_tax)

# or apply

tax_out <- lapply(rank_out$lineage[1:3], boldid_)
out_tbl <- data.frame(rank_out[1:3,], boldid = tax_out)

# then,

require(DT)

widget <- datatable(
  out_tbl, 
  escape=1,
  extensions = 'Buttons', options = list(
    pageLength = 25,  
    dom = 'Bfrtip',
    buttons = 
      list('copy', 'print', list(
        extend = 'collection',
        buttons = c('csv', 'excel'),
        text = 'Download'
      ))
    
  )
)

htmlwidgets::saveWidget(widget, paste0(tax.file, ".html"))

cat("\nDataTable conversion was done:", wang.taxonomy,"\n")

quit(save = "no")

