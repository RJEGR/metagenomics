# Identificar posicion de primers
# Esta estrategia es implementada por Cahalahan et al 
# https://benjjneb.github.io/dada2/ITS_workflow.html


rm(list = ls())

# test list :
path <- c('/Users/cigom/metagenomics/db/ncbi_nt_coi/CO1_COI_COX1_COXI_GENE_Eukaryota_nt/CO1_COI_COX1_COXI_GENE_Eukaryota_nr/pcr_seqs_screen/')

list_files <- list.files(path = path, pattern = '*.summary', full.names = TRUE)

str(Eukaryota_nr.pcr <- read.delim(list_files[1]))
str(Eukaryota_nr <- read.delim(list_files[2]))
str(Eukaryota_nr.unique <- read.delim(list_files[3]))


start <- Eukaryota_nr.pcr[, 'start']
end <- Eukaryota_nr.pcr[, 'end']
df <- as.data.frame(cbind(start,end))

poss <- vector()

i=1
for(i in 1:nrow(df)){
  poss <- append(poss, c(df[i,1]:df[i,2]))
  i+1
}

dens <- density(poss)

plot(dens)

# Viz blastn ----
require('ggplot2')

df <- dplyr::select(Eukaryota_nr.pcr, -seqname)

p <- ggplot() +
  geom_density(data = reshape2::melt(df),
               mapping = aes(x = value,
                             y = ..count.. 
               ), color = "grey") +
  labs(x = "Position", y = "Freq")

p + facet_wrap( ~ variable, scales = "free")


library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library("PrimerMiner")
packageVersion("PrimerMiner")


FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"  ## FORWARD &
REV <- "TANACYTCNGGRTGNCCRAARAAYCA"  ## REVERSE (inosina to N)
# REV <- "TAIACYTCIGGRTGICCRAARAAYCA"  ## REVERSE (inosina to N)

pathwd <- 'metagenomics/COI/run012/Ncyanomos/'
fasta.file <- paste0(pathwd, 'damisela_Ncyanomos_coi_short.fasta')

seqs <- readDNAStringSet(fasta.file)
seqs.width <- width(seqs)
names_seqs <- names(seqs)
names(seqs) <- seqs # This propagates to the tip labels of the tree


allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)



# 2

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFasta(fasta.file)), fixed = FALSE)
    return(sum(nhits > 0))
}
# Sanity check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fasta.file), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fasta.file), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fasta.file), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fasta.file))
# 

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor=NA)

out_align <- data.frame(alignment)

names(alignment) <- names_seqs

out_fasta <- c(rbind(paste0(">", names_seqs), out_align$alignment))

write(out_fasta)
file_name <- paste0(pathwd, "damisela_Ncyanomos_coi_short.align")
write(out_fasta, file=file_name)

# writePairwiseAlignments(alignment, file=paste0(pathwd, "damisela_Ncyanomos_coi_short.align"), Matrix=NA, block.width=50)
# clustering <- PrimerMiner::Clustering(fasta.file, vsearchpath = "Vsearch", id = 0.97, cmd = "" ,threshold = "Majority") # hay error en in file (con, 'r')
# Vsearch  -derep_fulllength damisela_Ncyanomos_coi_short.fasta -output /damisela_Ncyanomos_coi_short_drep.fasta


evaluate_primer()

evaluate_primer(alignment, FWD, forward = T, save = NULL, gap_NA = T, N_NA = T, mm_position
 = NULL, mm_type = NULL, adjacent = 2, sequ_names=T)

# 
# cd vsearch-2.12.0
# ./autogen.sh
# ./configure
# make
# make install 