dir <- '/Users/cigom/metagenomics/db/bold/BOLD_public_trim/'
tax <- read.delim(list.files(path = dir, pattern = 'bivalvia', full.names = T), sep = ' ', header = F) 
head(tax)

library(Biostrings)
ffile <- list.files(path = dir, pattern = 'BOLD_public_species.fasta', full.names = T)
dnaset <- Biostrings::readDNAStringSet(ffile)

dnaset <- dnaset[names(dnaset) %in% tax$V1]
dnaset <- DECIPHER::RemoveGaps(dnaset)
dnaset <- dnaset[width(dnaset) > 500]

# Dereplicate

dnaset <- sort(dnaset)
dnaset[!duplicated(dnaset)]

sum(tax[,1] %in% names(dnaset))
tax <- tax[tax[,1] %in% names(dnaset),]

tax$Class <- sapply(strsplit(tax[,2], ";"), `[`, 5)
tax$Order <- sapply(strsplit(tax[,2], ";"), `[`, 6)

table(tax$Order)

# convert to protein (decipher)
# and make tree (decipher) choosing genetic code of invertebrates
