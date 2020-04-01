
# Ricardo Gomez-Reyes, March, 2020
# 1. Read and combine data
# 2. Detect low-abundance (laf) features ----
# 3. Clean non-animal taxonomy ----
# 4. Filter laf and out file -----
# 5. Test visualizations

rm(list=ls())

# 0. Load functions ----
require(data.table)

url <- 'https://raw.githubusercontent.com/RJEGR/metagenomics/master/readtx.R'

source(url)

# set rank names ----


# la intension es tener un formato de rankins segun la base de datos.
# y hacer un check dbformat para en base de ello accesar al formato correcto.

dbformat <- c("silva", "bold", 'midori')
dbformat[dbformat %in% db]

frmt <- TRUE
dbformat = "bold"

# Silva ranking format
Ranks <- list(Reino=sym("Rank_1"), Filo=sym("Rank_2"), 
                Clase=sym("Rank_3"), Orden=sym("Rank_4"), 
                Familia=sym("Rank_5"), Genero=sym("Rank_6"),
                Especie=sym("Rank_7"))

# BOLD ranking format
Ranks <- list(root = sym("Rank_1"), Dominio=sym("Rank_2"), 
                Reino=sym("Rank_3"), Filo=sym("Rank_4"), 
             Clase=sym("Rank_5"), Orden=sym("Rank_6"), 
             Familia=sym("Rank_7"), Genero=sym("Rank_8"),
             Especie=sym("Rank_9"))

# Funcions (mover a ) ----
combine_features <- function(count_tbl, tax_tbl, k = 1, tax_header = TRUE) {
  
  require(data.table)
  require(dplyr)
  
  # 1. Read and combine data ----
  
  features <- data.table(read.table(count_tbl))
  tax <- data.table(read_rdp(tax_tbl, header = tax_header))
  size <- data.table(fsize = rowSums(features))
  nsams <- data.frame(nsams = c(apply(features, 1, function(x) {sum(x > 0)})))
  
  # test <- identical(rownames(features), rownames(tax))
  # 
  # if(test) { combine <- cbind(features, tax, size) }
  # else { combine <- features }
  
  combine <- cbind(features, tax, size, nsams)
  
  # 2. Detect low-abundance (laf) features ----
  
  combine[, ksize := ifelse(fsize > k,'Keep','Drop')]
  
  # 3.Clean non-animal taxonomy ----
  
  groups <- table(tax$Rank_1)
  
  g <- groups[groups == max(groups)]
  
  others <- combine[combine$Rank_1 != names(g),]
  others <- rename(others, !!!Ranks)
  
  combine <- combine[combine$Rank_1 == names(g),]
  
  nrow(combine) == g
  
  # 4. Filter laf and out file ---- 
  sam <- colnames(features)
  
  out <- filter(combine, ksize == 'Keep')
  
  out <- select(out, all_of(sam), !!!Ranks)
  
  
  # write.table(out, file = 'combine_tax_count_format.txt')
  
  return(list("clean"=out, "complete"=combine, "others" = others))
  
}


# SILVA ----
# Silva ranking format
Ranks <- list(Reino=sym("Rank_1"), Filo=sym("Rank_2"), 
              Clase=sym("Rank_3"), Orden=sym("Rank_4"), 
              Familia=sym("Rank_5"), Genero=sym("Rank_6"),
              Especie=sym("Rank_7"))

path <- "/Users/cigom/metagenomics/18S_MULTIRUN_XIXIMIS/downstream"
f <- list.files(pattern = "shared", full.names = TRUE, path = path)
t <- list.files(pattern = "taxonomy", full.names = TRUE, path = path)

out_silva <- combine_features(f, t, tax_header = TRUE)

write.table(out_silva$clean, file = paste0(path,'/','combine_tax_count_format1.txt' ))

# BOLD ----
path <- "/Users/cigom/metagenomics/COI/MULTIRUN/multirun_20190829_AMB_XIXIM/"

count_tbl <- list.files(pattern = "count.table", full.names = TRUE, path = path)

tax_tbl <- list.files(pattern = "taxonomy", full.names = TRUE, path = path)


out_bold <- combine_features(count_tbl, tax_tbl, tax_header = FALSE)

write.table(out_bold, file = paste0(path,'/','combine_tax_count_format2.txt' ))







# 5. Test visualizations -----
# para dar continuidad al combine, estaria bien meter algunas visualizaciones. etc.
others <- data.table(out_silva$others)
animalia <- data.table(out_silva$complete)

table(others$Resolution)

aglom_ab <- function(x, rank) {
  x <- as.data.frame(x)
  tax_g <- aggregate(x[,'fsize'],  by = list(x[, rank]), FUN = sum)
  n_asvs <- aggregate(x[, rank],  by = list(x[, rank]), FUN = length)
  
  if(identical(tax_g[,1], n_asvs[,1])){
      tax_out <- data.frame(lineage = tax_g[,1], 
                            Abundance = tax_g[,2], 
                            Size = n_asvs[,2])
  }
  return(tax_out)
}

aglom_ab(others, rank = 'Genero')

# obtener un indice de 'peso' que nos deje saber la contribucion de los features para alcanzar dicha abundancia. pj fsize y nsams de estos 7 features del mismo genero

others[others$Genero == 'Acanthometra']

