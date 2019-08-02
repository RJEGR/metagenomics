# 
# Agosto 2019
# Ricardo Gomez

rm(list = ls())

# Functions ----
dada2 <- function(x, OMEGA_A, err) {
  dada(x, err = err, 
       multithread = threads, 
       selfConsist = FALSE, OMEGA_A = OMEGA_A)
}

learn_Err <- function(x) { learnErrors(x, 
                                       nbases = 1e+08,
                                       MAX_CONSIST = MAX_CONSIST,
                                       multithread = threads, 
                                       OMEGA_C = 0, 
                                       errorEstimationFunction = loessErrfun,
                                       verbose = TRUE,
                                       randomize = FALSE)
}



# setdir ----
dir <- '/Users/cigom/metagenomics/COI/run15'
setwd(dir)
# 
## Checking and Load packages ----
# 
.cran_packages <- c('dplyr', 'purrr', 'tibble', 'reshape2', 'ggplot2', 'tidyr')
.bioc_packages <- c("dada2")


.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)


packageVersion("dada2")

## Variables to check/change ----
#

# Do trimming
do_at_trim <- TRUE
# Do Filtering
do_derep <- TRUE
# Do Assign Taxonomy
do_at_assig <- FALSE

source("config.R")


# OMEGA_C <- 0 # which means that all reads are counted and contribute to estimating the error rates
threads = TRUE

# OMEGA_A <- 1e-100
# KDIST_CUTOFF <- 0.42
# setDadaOpt(OMEGA_A = as.double(OMEGA_A)) # KDIST_CUTOFF = KDIST_CUTOFF


data_path = dir
out_path = paste0(dir, '/dada2_asv')


## DADA2 pipeline ---- MOCK_HITS
fnFs <- sort(list.files(data_path, pattern = "015-Mock27-COI-Zoo_S56_L001_R1_001.fastq.gz", full.names = T))
fnRs <- sort(list.files(data_path, pattern = "015-Mock27-COI-Zoo_S56_L001_R2_001.fastq.gz", full.names = T))

plotQP(list.files(data_path, pattern = 'fastq.gz', full.names = TRUE))

filtFs <- file.path(out_path, paste0("015-Mock27-COI-Zoo_F_filt.fastq"))
filtRs <- file.path(out_path, paste0("015-Mock27-COI-Zoo_R_filt.fastq"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = fat_trunclen, 
                     truncQ = fat_trunQ, 
                     maxEE= fat_maxee, 
                     minLen = fat_minlen,
                     trimLeft = 26, maxN = 0, 
                     rm.phix = rm_phix, compress = T, multithread = T)

plotQP(list.files(out_path, pattern = 'filt.fastq', full.names = TRUE))

# Dereplication ----
derepFs <- derepFastq(filtFs, verbose = T)
derepRs <- derepFastq(filtRs, verbose = T)

# Learn Error ----
MAX_CONSIST <- 20

errF <- learn_Err(filtFs)
errR <- learn_Err(filtRs)

# plotErrors(errF, nominalQ=TRUE)

# dada algorithm
OMEGA_A <- 1e-120

dadaFs <- dada2(derepFs, OMEGA_A, errF)
dadaRs <- dada2(derepRs, OMEGA_A, errF)

# merging ----
propagateCol <- names(dadaFs)

nrow(mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                           minOverlap = 12, maxMismatch = 1,
                           verbose=TRUE, returnRejects = TRUE, 
                           propagateCol=propagateCol))



dim(seqtab <- makeSequenceTable(mergers))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = T)

# write data ----
asv_seqs <- colnames(seqtab.nochim)
n_asv <- dim(seqtab.nochim)[2]
asv_headers <- vector(n_asv, mode="character")

for (i in 1:n_asv) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}


asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
colnames(asv_tab) <- 'Mock'

out_prefix <- paste0('OMEGA_A_', gsub("-", "_", OMEGA_A), '_')

write.table(asv_tab,
            file=paste0(out_path, "/", 
                        out_prefix, "ASVs_count.table"), 
            sep="\t", 
            row.names = TRUE, 
            col.names = TRUE,
            quote=FALSE
)

# and fasta
save_fasta <- c(rbind(asv_headers, asv_seqs))
write(save_fasta, file=paste0(out_path, "/", out_prefix, "asvs.fasta"))

# quit(save = 'no')

# Taxonomy summ ----
wang.taxonomy <- 'OMEGA_A_1e_120_asvs.ALL.wang.taxonomy'
count_tbl <- 'OMEGA_A_1e_120_ASVs_count.table'
fasta_file <- 'OMEGA_A_1e_120_asvs.fasta'
path <- '/Users/cigom/metagenomics/COI/run15/dada2_asv'

# Load files
source(file = "~/Documents/GitHub/metagenomics/readtx.R")
tax.file <- list.files(path = path, full.names = TRUE, pattern = wang.taxonomy)[1]
count.file <- list.files(path = path, full.names = TRUE, pattern = count_tbl)
fasta.file <- list.files(path = path, full.names = TRUE, pattern = fasta_file)

TL <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

scale2 <- c("Domain"="#edf8b1",  "Kingdom"="#7fcdbb", "Phylum"="#2c7fb8",  
            "Class"="#feb24c",  "Order"="#addd8e",  "Family"="#31a354",
            "Genus"="#bcbddc", "Species"="#756bb1")

tax <- read_rdp(tax.file)

colnames(tax) <- c(TL, 'SL')

tax <- data.frame(ASV = rownames(tax), tax)

out0 <- bbold_(tax, fasta.file, count.file,  rel_ab = FALSE)

out0[is.na(out0)] <- 'Undetermined_R'

# data.frame(table(out0$Species))
# Alluvial ----
alluv_in <- out0[,TL[4:9]]

library(ggalluvial)
alluv_long <- to_lodes_form(data.frame(alluv_in), key = "Rank",axes = 1:6)
alluv_long <- filter(alluv_long, stratum != 'Undetermined_R')


is_alluvia_form(alluv_long, silent = TRUE)

ggplot(data = alluv_long,
       aes(x = Rank, stratum = stratum, alluvium = alluvium,
           label = stratum)) +
  geom_stratum() + 
  geom_text(stat = "stratum", size = 2.5) +
  geom_flow( stat = "alluvium", # aes(fill = alluvium),
            aes.bind = TRUE, lode.guidance = "rightward") +
  theme_minimal() +
  ggtitle("The frequency distribution of zooplankton in the Mock-015") +
  xlab("Level of resolution") + ylab("Number of ASVs")

#
agg_species <- aglom_ab(out0, 'Species')

