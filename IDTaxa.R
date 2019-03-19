
# Author: Ricardo Gomez-Reyes | modified from: Erik S. Wright (2018) - Classify Sequences

# The classification process is split into two parts: 
# training carried out by:
# LearnTaxa and testing with IdTaxa
# How to use: Rscript --vanilla IDTaxa.R input[string] reference[string] taxonomy[string] threshold[interger - by default 99]
# reference are fasta format of sequence references 
# taxonomy are assignments for sequences embed in the reference file; two columns as in rdp-mothur format (';' separated ranks-by-rank)
## >> Note: All elements of taxonomy must contain 'Root;'

# Rscript --vanilla IDTaxa.R run012_relax_ASVs.P68.worms.fasta BOLD_public_species.fasta BOLD_public_species.tax.Root 80
## Clean workspace
rm(list=ls()); 

# Path ----
library(DECIPHER)

path = '/Users/cigom/metagenomics/COI/run012'
setwd(path)

##  args insertion ----

args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
  stop("!!!\n. 
       PLEASE, INPUT NECESSARY FILES IN THE SYNTAXIS AS FOLLOW EXAMPLE:\n
       Rscript --vanilla IDTaxa.R input[string] reference[string] taxonomy[string] threshold[interger] .\n", call.=FALSE)
} else { 
  query <- args[1]
  seqs_name <- args[2]
  tax_name <- args[3] }

# Read inputs ----

# query <- "run012_relax_ASVs.P68.worms.fasta"
q_path <- paste0(path, "/", query)

# Sequnces references path
# seqs_name <- "IctiosConsenso.fasta"
# seqs_name <- 'BOLD_public_species.fasta'
seqs_path <- paste0(path, "/", seqs_name)

# taxonomy assigments path
# tax_name <- 'IctiosConsenso.tax'
# tax_name <- 'BOLD_public_species.tax.Root'
tax_path <- paste0(path, "/", tax_name)

#
tag <- strsplit(query, "[.]")[[1]]
tag_ref <- strsplit(seqs_name, "[.]")[[1]]
out_name <- paste0(paste(tag[-length(tag)], collapse="."), "_",
                   paste(tag_ref[-length(tag_ref)], collapse="."), "_",
                                                  "DECIPHER.taxonomy")
#
# Parameters ----
# 
if (length(args) > 3) { threshold <- as.numeric(args[4]) } else { threshold <- 99 }


## Training parameters

maxGroupSize <- 1 # max sequences per label (>= 1); can be set to Inf (infinity) to allow for an unlimited number of sequences per group.
maxIterations <- 1 # must be >= 1; also can remove mislabeled in the training data by turning itinerations > 1 
allowGroupRemoval <- FALSE # if maxIterations > 1 set allowGroupRemoval <- TRUE


#
# Load sequence from reference data-base: ----
#
# Read sequences into memory
seqs <- readDNAStringSet(seqs_path) # or  readRNAStringSet for RNA
# if gaps in sequences:
seqs <- RemoveGaps(seqs)

#
# Read the taxonomic assigments from data-base ----
# 
taxonomy <- read.table(tax_path, header = FALSE, stringsAsFactors = FALSE)

# Renombramos el nombre de secuencias con la asignacion correspondiente:
if (identical(taxonomy$V1, names(seqs))) {
               names(seqs) <- taxonomy$V2
               groups <- names(seqs)
               groups <- gsub("(.*)(Root;)", "\\2", groups)
               groupCounts <- table(groups)
               cat("number of groups in reference are:\n", length(u_groups <- names(groupCounts)))
    } else {
               groups <- taxonomy$V2
               groups <- gsub("(.*)(Root;)", "\\2", groups)
               groupCounts <- table(groups)
               cat("number of groups in reference are:\n", length(u_groups <- names(groupCounts)))
               }
                 
#
# Process the training set of data-base: ----
#

# Pruning the traning set
# Count the number of representative per group

# maxGroupSize <- 2 # max sequences per label (>= 1); can be set to Inf (infinity) to allow for an unlimited number of sequences per group.
remove <- logical(length(seqs))

for (i in which(groupCounts > maxGroupSize)) {
  index <- which(groups==u_groups[i])
  keep <- sample(length(index),
                 maxGroupSize)
  remove[index[-keep]] <- TRUE
}

sum(remove) # number of sequences eliminated

# Iteratively training the classifier
# Is it will identify any training sequences whose 
# assigned classifications completely (with very high confidence) disagree with their predicted classification.

# maxIterations <- 1 # must be >= 1; also can remove mislabeled in the training data by turning itinerations > 1 
# allowGroupRemoval <- FALSE

probSeqsPrev <- integer()
taxid <- NULL

for (i in seq_len(maxIterations)) {
  cat("Training iteration: ", i, "\n", sep="")
  # train the classifier
  trainingSet <- LearnTaxa(seqs[!remove],
                           names(seqs)[!remove],
                           taxid)
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  if (i==maxIterations)
    break
  probSeqsPrev <- probSeqs
  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
}


# View training results:

trainingSet
plot(trainingSet)

#
# Classifying Sequences: ----
# 

q_seqs <- readDNAStringSet(q_path)
q_seqs <- RemoveGaps(q_seqs) # if gaps remove it

ids <- IdTaxa(q_seqs, 
              trainingSet, 
              type="extended",
              strand = "top",
              threshold = threshold,
              processors = 1)

#ids[1:5]
#ids[[2]]
#ids[c(10, 25)]
#c(ids[10], ids[25])

plot(ids, trainingSet)

# subset any rank
rank <- sapply(ids,
                 function(x) {
                   w <- which(x$rank=="Root")
                   if (length(w) != 1) {
                     "unknown"
                   } else {
                     x$taxon[w]
                   }
                 })
# table(rank)

taxon <- sapply(ids,
                  function(x)
                    x$taxon[length(x$taxon)])

# Output assignments: ----
output <- data.frame(sapply(ids,
                     function(x)
                       paste(x$taxon,
                             collapse=";")))


output <- data.frame(ASV = rownames(output), Taxonomy = output[,1])
write.table(output, paste0(path, '/', out_name), sep = "\t", row.names = FALSE, col.names = TRUE)
      

#
output_confidence <- sapply(ids,
                 function (id) {
                   paste(id$taxon,
                         " (",
                         round(id$confidence, digits=1),
                         "%)",
                         sep="",
                         collapse="; ")})


writeLines(output_confidence, paste0(path, '/', out_name, ".confidence"))


quit(save='no')

summary(confidence <- sapply(ids,
                             function (id) {
                               round(id$confidence, digits=1)
                             }))
