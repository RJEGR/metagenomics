
.cran_packages <- c('dplyr', 'purrr', 'tibble', 'reshape2')
.bioc_packages <- c("metacoder")

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


# =======
# INPUTS:
# =======
# sample files
print(hmp_samples) # SAMPLE FACTOR FILE
print(hmp_otus) # COUNT_TABLE W/ LINEAGE AND ASV ID

# One challenge this data presents is the format of the taxonomic information.
hmp_otus$lineage[1:4]

obj <- parse_tax_data(hmp_otus, # Input all the count_table curated with ASV_ID and lineage assignation
                      class_cols = "lineage", # the column that contains taxonomic information
                      class_sep = ";", # The character used to separate taxa in the classification
                      class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                      class_key = c(tax_rank = "info", # A key describing each regex capture group
                                    tax_name = "taxon_name"))

#

