
# Ricardo Gomez-Reyes, March, 2020
# 1. Read and combine data
# 2. Detect low-abundance (laf) features ----
# 3. Clean non-animal taxonomy ----
# 4. Filter laf and out file -----
# 5. Test visualizations

# ==============
## Checking and Load packages ----
# ==============

.cran_packages <- c('tidyverse', 'purrr','ggplot2', 'RColorBrewer')

.bioc_packages <- c("Biostrings", "IRanges")

# 1.
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
# 2.
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

rm(list=ls())

# 0. Load functions ----
url <- 'https://raw.githubusercontent.com/RJEGR/metagenomics/master/readtx.R'
source(url)

path <- "/Users/cigom/metagenomics/Franccesco/callahan_multirun/outputs/"

# pathnames ----
# f_name <- "dada_ASVs_count.table"
# t_name <- "dada_ASVs.w2pr2_worms_API02.wang.taxonomy"
# seq_name <- 'dada_ASVs.fasta'

f_name <- "dada_fra_ASVs_count.table"
t_name <- "dada_fra_ASVs.w2pr2_worms_API02.wang.taxonomy"
seq_name <- 'dada_fra_ASVs.fasta'

fasta_file <- list.files(pattern = seq_name, full.names = TRUE, path = path)
count_tbl <- list.files(pattern = f_name, full.names = TRUE, path = path)
tax_tbl <- list.files(pattern = t_name, full.names = TRUE, path = path)


# set rank names ----


# la intension es tener un formato de rankins segun la base de datos.
# y hacer un check dbformat para en base de ello accesar al formato correcto.
db <- 
dbformat <- c("silva", "bold", 'midori')
dbformat[dbformat %in% db]

frmt <- TRUE
dbformat = "midori"
# SILVA ----
# Silva ranking format

# Silva ranking format
Ranks <- list(Reino=sym("Rank_1"), Filo=sym("Rank_2"), 
                Clase=sym("Rank_3"), Orden=sym("Rank_4"), 
                Familia=sym("Rank_5"), Genero=sym("Rank_6"),
                Especie=sym("Rank_7"))

# # BOLD ranking format
# Ranks <- list(root = sym("Rank_1"), Dominio=sym("Rank_2"), 
#                 Reino=sym("Rank_3"), Filo=sym("Rank_4"), 
#              Clase=sym("Rank_5"), Orden=sym("Rank_6"), 
#              Familia=sym("Rank_7"), Genero=sym("Rank_8"),
#              Especie=sym("Rank_9"))

# Funcions (mover a ) ----
combine_features <- function(count_tbl, tax_tbl, 
                             fasta_file, k = 1, tax_header = F) {
  
  url <- 'https://raw.githubusercontent.com/RJEGR/metagenomics/master/readtx.R'
  
  source(url)
  require(dplyr)
  
  # 1. Read and combine data ----
  
  features <- read.table(count_tbl, row.names = 1)
  tax <- read_rdp(tax_tbl, header = tax_header)
  seqs <- readDNAStringSet(fasta_file)
  
  tax_resolution <- tax$Resolution
  
  # select tax lineage ----
  taxRrank <- colnames(tax)
  
  tax <- tax[taxRrank %in% Ranks]
  
  taxRrank <- names(Ranks[taxRrank %in% Ranks])
  colnames(tax) <- taxRrank
  
  # tax <- mutate_all(tax, funs(str_replace_all(., c("Undetermined_R"=NA, " "="_"))))
  
  size <- tibble(fsize = rowSums(features))
  nsams <- data.frame(nsams = c(apply(features, 1, function(x) {sum(x > 0)})))
  
  datavis <- data.frame(fsize = rowSums(features),
                  nsams = c(apply(features, 1, function(x) {sum(x > 0)})),
                  res = tax_resolution, Rank1 = tax[,1])
  
  pct <- function(x) {x / sum(x) * 100}
  
  datavis %>%
    mutate(wrap = ifelse(log10(fsize) > 2, 'a', 'b'),
           pct = pct(fsize)) %>%
    ggplot(aes(nsams, log10(fsize), 
               color = as.factor(Rank1),
               size = pct)) +
    geom_point(alpha = 0.7, shape = 16) +
    theme(legend.position = "top") +
    scale_color_brewer(type='qual', palette=2) +
    theme_classic() +
    facet_grid(wrap ~., scales = 'free_y', space="free") + 
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    labs(color = 'Resolution', 
         x = 'Frequency in samples', y = 'log(Number of sequences)')
  
  test <- identical(rownames(features), rownames(tax))

  if(test) { combine <- cbind(features, tax) } else { combine <- features }
  
  rownames(combine) <- rownames(features)
  
  #combine <- cbind(features, tax, size, nsams)
  
  # 2. Detect low-abundance (laf) features ----
  combine %>%
    as_tibble(rownames = 'ids') %>%
    mutate(fsize = rowSums(features)) %>%
    mutate(ksize = ifelse(fsize > k,'Keep','Drop')) -> combine
  
  # 3.Clean non-animal taxonomy ----
  
  groups <- table(tax[,1])
  
  g <- groups[groups == max(groups)]
  
  others <- combine[combine[, taxRrank[1]] != names(g),]
  
  combine <- combine[combine[, taxRrank[1]] == names(g),]
  
  nrow(combine) == g
  
  # 4. Filter laf and out file ---- 
  sam <- colnames(features)

  
  out <- filter(combine, ksize == 'Keep')
  
  out %>%
    select_at(vars(ids,sam)) -> features_filtered
  
  out %>%
    select_at(vars(all_of(taxRrank))) %>%
    unite(sep = ";", remove = TRUE, col = 'Taxonomy') %>%
    cbind(out$ids, .) -> tax_filtered
  
  ids_filtered <- out$ids 
  
  seqs <- seqs[names(seqs) %in% ids_filtered,]
  
  # write.table(out, file = 'combine_tax_count_format.txt')
  
  
  return(list("ff" = features_filtered, "tf" = tax_filtered, 
              "others" = others, "seqs" = seqs))
  
}


out <- combine_features(count_tbl, tax_tbl, fasta_file, tax_header = FALSE)

system(paste0('mkdir -p ', path, '/filtered_fra'))

outpath <- paste0(path, "/filtered_fra")

write.table(out$ff, 
            file = paste0(outpath,'/',  basename(count_tbl)), 
            row.names = F, sep = '\t', quote = F)

write.table(out$tf, 
            file = paste0(outpath,'/',  basename(tax_tbl)), 
            row.names = F, sep = '\t', col.names = F, quote = F)

writeXStringSet(out$seqs, 
                paste0(outpath,'/',  basename(fasta_file)))

ggsave(filename = "log_plot.png", dpi = 200 ,path = outpath)






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

