# Identificar posicion de primers
# Esta estrategia es implementada por Cahalahan et al 
# https://benjjneb.github.io/dada2/ITS_workflow.html

library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library("PrimerMiner")
packageVersion("PrimerMiner")


FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"  ## FORWARD &
REV <- "TAIACYTCIGGRTGICCRAARAAYCA"  ## REVERSE (inosina to N)

fasta.file <- 'damisela_Ncyanomos_coi_short.fasta'

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

primerHits <- function(primer, seqs) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFasta(fasta.file)), fixed = FALSE)
    return(sum(nhits > 0))
}
# Sanity check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, seqs = seqs), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, seqs = seqs), 
    REV.ForwardReads = sapply(REV.orients, primerHits, seqs = seqs), 
    REV.ReverseReads = sapply(REV.orients, primerHits, seqs = seqs))

# 
seqs <- readDNAStringSet(fasta.file)
seqs.width <- width(seqs)
names_seqs <- names(seqs)
names(seqs) <- seqs # This propagates to the tip labels of the tree

alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), anchor=NA)

asv_fasta <- c(rbind(asv_headers, asv_seqs))

writePairwiseAlignments(alignment, file="damisela_Ncyanomos_coi_short.align", Matrix=NA, block.width=50)
# clustering <- PrimerMiner::Clustering(fasta.file, vsearchpath = "Vsearch", id = 0.97, cmd = "" ,threshold = "Majority") # hay error en in file (con, 'r')
# Vsearch  -derep_fulllength damisela_Ncyanomos_coi_short.fasta -output /damisela_Ncyanomos_coi_short_drep.fasta


evaluate_primer()




evaluate_primer(alignment, FWD, forward = T, save = NULL, gap_NA = T, N_NA = T, mm_position
 = NULL, mm_type = NULL, adjacent = 2, sequ_names=T)


cd vsearch-2.12.0
./autogen.sh
./configure
make
make install 