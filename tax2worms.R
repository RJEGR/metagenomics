
# 
# Junio 2019
# Ricardo Gomez
# Rscript --vanilla ~/Documents/GitHub/metagenomics/tax2html.R wang.taxonomy count_tbl fasta_file

rm(list = ls())

options(stringsAsFactors = FALSE)

wang.taxonomy <-  'taxonomy_80_COI_filtered.csv' # 'multirun_ASVs.ALL.wang.taxonomy'
count_tbl <- 'feature_table_80_COI_filtered.csv'# 'multirun_ASVs_count.table'

fasta_file <- 'multirun_ASVs.fasta'

rank <- 'Species'

# ------ functions

# Include the convertion to retrive BOLD id or BOLD BIND with

boldid_ <- function(x) {
  require(taxize)
  x <- unique(x)
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

boldid_3 <- function(x) {
  x_ <- NULL
  for (i in 1:nrow(x)) {
    x_[[i]] <- boldid_(x$lineage[i])
  }
  return(x_)
}

names2wormsid <- function(mynames, accepted = FALSE) {
  worms <- get_wormsid(mynames, accepted = accepted, ask = FALSE)
  
  id <- worms[1:length(worms)]
  uri <- attributes(worms)$uri
  
  out <- data.frame(linage = mynames, wormsid = id, uri = uri )
  
  return(out)
}

# wormsid2name <- function(){}

saveWget <- function(out_tbl, pageLength) {
  require(DT)
  
  widget <- datatable(
    out_tbl, 
    escape=1,
    extensions = 'Buttons', options = list(
      pageLength = pageLength,  
      dom = 'Bfrtip',
      buttons = 
        list('copy', 'print', list(
          extend = 'collection',
          buttons = c('csv', 'excel'),
          text = 'Download'
        ))
      
    )
  )
  
  return(widget)
  
}

# Path ----
# path <- getwd()

path <- '/Users/cigom/metagenomics/Franccesco/multirun_xixim_coi/'

# Load files ----

tax.file <- list.files(path = path, full.names = TRUE, pattern = wang.taxonomy)
count.file <- list.files(path = path, full.names = TRUE, pattern = count_tbl)
fasta.file <- list.files(path = path, full.names = TRUE, pattern = fasta_file)


cat('File to process: ', tax.file[1], "\n")
cat('with ', count_tbl[1], 'and ', fasta_file)

TL <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ----- 

source(file = "~/Documents/GitHub/metagenomics/readtx.R")

tax <- read_rdp(tax.file[1])
colnames(tax) <- c(TL, 'SL')


cat('Sample size is:', nrow(tax), "\n")

# la funcion boldid_2 toma el ultimo rank asignado por asv y hace una busqueda del boldid (tiempo demandante)

# Esta version es mas rapida (boldid_3), resume el linaje al nivel de interes, en este ejemplo se usa especie 
# ademas se incluye la abundacia aglomerada (sea numero de reads o abundancia relativa)

tax <- data.frame(ASV = rownames(tax), tax)

out0 <- bbold_(tax, fasta.file[1], count.file[1],  rel_ab = FALSE)

str(rank_out <- data.frame(aglom_ab(out0, rank)))

# connect to a secured network (or error with Web authent)

bold3_tax <- boldid_3(rank_out) # time demand, make uniques the ranks!
out_tbl <- data.frame(rank_out, boldid = bold3_tax)

widget <- saveWget(out_tbl, 20)

htmlwidgets::saveWidget(widget, paste0(tax.file[1], "." ,rank, ".html"))

cat("\nDataTable conversion was done for:", wang.taxonomy,"\n", "in output: ", paste0(tax.file[1], ".html"))


u_names <- sort(out_tbl$lineage)

# 1 non-interactive worms name selection 



# children(worms$linage, db = 'worms')

library(worrms)

end_date <- format(Sys.Date() - 1, "%Y-%m-%dT%H:%M:%S+00:00")
# end_date_bkp <- "2019-11-07T00:00:00+00:00"
start_date <- format(Sys.Date() - 365*10, "%Y-%m-%dT%H:%M:%S+00:00")

wm_records <- wm_records_date(start_date = start_date, end_date = end_date, marine_only = TRUE)

# wm_records_names(worms$linage, marine_only = TRUE)

# OR 

wm_records <- read.table('~/Documents/GitHub/metagenomics/converted2worms/midori2worms_conversion.tab.2019Mar06')

# generar la base de records de worms como arriba y recuperar los ids de worms correspondientes a nuestros resultados!
# ref https://docs.ropensci.org/worrms/
worms
wm_records
inner_join()

# quit(save = "no")
# boldid_(tax$Phylum)

tax2name <- select(tax, names(tax)[5:10])

query <- unique(tax2name$Phylum)

# HAY UN PROBLEMON PARA PROCESAR DATOS , SE REPITEN LOS NOMBRES EN EL MATCHING DE NOMBRES, EJ Nematoda ES ENCONTRADO COMO Nematoda Y 492652-

# wid <- get_wormsid_(query, searchtype = "scientific", accepted = TRUE, ask = TRUE)
# wid <- do.call(rbind, wid)
# id <- wid[, 'AphiaID']

#gid <- get_wormsid(query, searchtype = "scientific")
#match <- attributes(gid)$match
# id <- gid[1]

out <- names2worms_(query)




tax2name <- boldid_2(tax2name)

worms <- names2wormsid(head(u_names), accepted = TRUE)


