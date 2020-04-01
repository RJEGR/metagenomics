rm(list=ls())

# set paths and functions ----

path <- c("/Users/cigom/metagenomics/MG_18S/run017_M4/")

url <- 'https://raw.githubusercontent.com/RJEGR/metagenomics/master/readtx.R'

source(url)

Ranks <- c("Reino", "Filo", "Clase", "Orden", "Familia", "Genero", "Especie")

makeMT <- function(dat) {
  samples <- names(dat)
  run <- sapply(strsplit(samples, "_"), `[`, 1)
  cruise <- sapply(strsplit(samples, "_"), `[`, 2)
  subject <- sapply(strsplit(samples, "_"), `[`, 3)
  marker <- sapply(strsplit(samples, "_"), `[`, 4)
  sample_type <- sapply(strsplit(samples, "_"), `[`, 5)
  transect <- substr(subject,1,1)
  
  out <- data.frame(Corrida = run, 
                    Crucero=cruise, 
                    Estación=subject, 
                    Transecto=transect, 
                    muestra=sample_type, 
                    marcador =  marker)
  
  rownames(out) <- samples
  return(out)
  
}


obj.dir <- path
setwd(obj.dir)

# # # # # #
# load data (after processing) ----

mcount <- list.files(obj.dir, full.names = TRUE, pattern = "*.shared")
mtax <- list.files(obj.dir, full.names = TRUE, pattern = "*.taxonomy")

dat <- read.table(mcount, header = T)
tax <- read_rdp(mtax, header = T)

# for some reason split columns to rows in count_tbl
# headers: label, Group, numOtus, Otu00001 .... Otu0000N

samples <- dat$Group
dat <- data.frame(t(dat[-c(1:3)]))
names(dat) <- samples

# Load metadata

samples <- makeMT(dat)

# label tax lineage
names(tax) <- Ranks
tax <- tax[Ranks]


# construct phyloseq object to clean data

identical(names(dat), rownames(samples))
identical(rownames(dat), rownames(tax))

phyloseq = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')), 
                    sample_data(samples))


# 

physeq <- subset_taxa(phyloseq, Reino == "Animalia")
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)

#keepTaxa = apply(X = as(otu_table(physeq), "matrix") > 0, 
#                 MARGIN = 1, FUN = sum) > 2 

#table(keepTaxa)
# clean.phyloseq = prune_taxa(keepTaxa, physeq)


write.table(otu_table(physeq), 
            paste0(mcount, '_animalia_feature_count.table'))

# Rscript --vanilla estimate_richness.R *_animalia_feature_count.table

quit(save = 'no')

