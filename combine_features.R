
# Ricardo Gomez-Reyes

# 0. Load functions ----
require(data.table)

path <- "/Users/cigom/metagenomics/18S_MULTIRUN_XIXIMIS/downstream"

setwd(path)

url <- 'https://raw.githubusercontent.com/RJEGR/metagenomics/master/readtx.R'

source(url)
# rank ----
silva <- TRUE

if(silva) {
  # Silva Ranks
  Ranks <- c("Reino", "Filo", "Clase", "Orden", 
              "Familia", "Genero", "Especie")
} else {
  # BOLD ranking 
  Ranks <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
}

# 1. Read and combine data ----

f <- list.files(pattern = "shared", full.names = TRUE)
t <- list.files(pattern = "taxonomy", full.names = TRUE)

features <- data.table(read.table(f))
tax <- data.table(read_rdp(t, header = TRUE))

size <- data.table(fsize = rowSums(features))
test <- identical(rownames(features), rownames(tax))


if(test) { combine <- cbind(features, tax, size) } else { combine <- features }

# 2.Clean non-animal taxonomy ----

groups <- table(tax$Rank_1)

g <- groups[groups == max(groups)]

others <- combine[combine$Rank_1 != names(g),]

combine <- combine[combine$Rank_1 == names(g),]

nrow(combine) == g

# 3. Check low-abundance (laf) features ----

k <- 1

combine[, ksize := ifelse(fsize > k,'Keep','Drop')]

# x. agglomerate ??

# 4. Filter laf and out file 

sam <- c(colnames(features), Ranks)

out <- filter(combine, ksize == 'Keep')
out <- select(combine, all_of(sam))

write.table(out, file = 'combine_tax_count_format.txt')
