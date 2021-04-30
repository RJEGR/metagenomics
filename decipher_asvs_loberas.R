rm(list = ls())

.bioc_packages <- c("Biostrings","DECIPHER") # "phangorn"

# Load packages into session, and print package version
sapply(.bioc_packages, require, character.only = TRUE)

# file <- dir(pattern="centroids_ab.fa")[1]

ssu_path <- '/Users/cigom/metagenomics/Loberas_MG/fastq/dada2_asv/ssu_pick_dir'
fastaN <- 'ssu_pick_dir.bacteria.fa'

file <- list.files(path = ssu_path, pattern = fastaN, full.names = T)
Reffile <- '/Users/cigom/metagenomics/Loberas_MG/Streptococcus mutans_ref_Seq.fasta'

dna <- readDNAStringSet(file, format="fasta")

dnaRef <- readDNAStringSet(Reffile)

dna <- c(dnaRef, dna)

aligment <- AlignSeqs(dna)
BrowseSeqs(aligment)

# name "Vertebrate Mitochondrial" ,
GCode <- getGeneticCode("SGC0")

gT <- lapply(order(width(dna), decreasing=TRUE),
             function(x) {
               attr(x, "height") <- 0
               attr(x, "label") <- names(dna)[x]
               attr(x, "members") <- 1L
               attr(x, "leaf") <- TRUE
               x
             })


attr(gT, "height") <- 0.5
attr(gT, "members") <- length(dna)
class(gT) <- "dendrogram"

# use the guide tree as input for alignment


DNA <- AlignTranslation(dna,
                        guideTree=gT,
                        iterations=0,
                        refinements=0,
                        geneticCode = GCode)

DNA <- AdjustAlignment(DNA)

BrowseSeqs(DNA)

# 

DNA <- DNA[!names(DNA) %in% cleanSeqs]

writeXStringSet(DNA, filepath = paste0(file, ".decipher.afa"), append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")


