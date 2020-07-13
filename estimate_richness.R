#!/usr/bin/env Rscript

# Use:  
# Rscript --vanilla estimate_richness.R ASVs_count.table measures
# measures:
# Observed-Chao1-ACE-Fisher-Shannon-InvSimpson

# pks required
# R version 3.5.0
# vegan - 2.5.6
# knitr - 1.25

# ================
# Defining outputs
# ================
outpath <- getwd()

args = commandArgs(trailingOnly=TRUE)

file = args[1]

if (is.na(args[2])) {
  measures <- NULL
} else 
  measures = unlist(lapply(strsplit(args[2], "-", fixed=TRUE), 
                                          function(x) return(x)))


cat("Using ", file, " as input file \n And use: ",
       measures, "\nto calculate diversity")

estimate_richness <- function (feature_tbl, split = TRUE, measures = NULL) 
{
  if (!any(feature_tbl == 1)) {
    warning("The data you have provided does not have\n", 
            "any singletons. This is highly suspicious. Results of richness\n", 
            "estimates (for example) are probably unreliable, or wrong, if you have already\n", 
            "trimmed low-abundance taxa from the data.\n", "\n", 
            "We recommended that you find the un-trimmed data and retry.")
  }
  if (!split) {
    OTU <- rowSums(feature_tbl)
  }
  else if (split) {
    OTU <- as(feature_tbl, "matrix")
      OTU <- t(OTU)
  }
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", 
                "InvSimpson", "Fisher")
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", 
                        "simpson", "invsimpson", "fisher")
  if (is.null(measures)) {
    measures = as.character(renamevec)
  }
  if (any(measures %in% names(renamevec))) {
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% 
                                                            measures]
  }
  if (!any(measures %in% renamevec)) {
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  outlist = vector("list")
  estimRmeas = c("Chao1", "Observed", "ACE")
  if (any(estimRmeas %in% measures)) {
    require(vegan)
    outlist <- c(outlist, list(t(data.frame(estimateR(OTU)))))
  }
  if ("Shannon" %in% measures) {
    outlist <- c(outlist, list(shannon = diversity(OTU, 
                                                   index = "shannon")))
  }
  if ("Simpson" %in% measures) {
    outlist <- c(outlist, list(simpson = diversity(OTU, 
                                                   index = "simpson")))
  }
  if ("InvSimpson" %in% measures) {
    outlist <- c(outlist, list(invsimpson = diversity(OTU, 
                                                      index = "invsimpson")))
  }
  if ("Fisher" %in% measures) {
    fisher = tryCatch(fisher.alpha(OTU, se = TRUE), warning = function(w) {
      warning("Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
      suppressWarnings(fisher.alpha(OTU, se = TRUE)[, 
                                                    c("alpha", "se")])
    })
    if (!is.null(dim(fisher))) {
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    }
    else {
      outlist <- c(outlist, Fisher = list(fisher))
    }
  }
  out = do.call("cbind", outlist)
  namechange = intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, 
                   colnames(out), ignore.case = TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop = FALSE]
  out <- as.data.frame(out)
  return(out)
}

# 

require(dplyr)

feature_tbl <- read.table(file, header = T) %>%
  select_if(., is.numeric)

out <- estimate_richness(feature_tbl, measures = measures)

write.table(out, file = paste0(outpath, "/estimated_richness.csv"))

knitr::kable(out, caption = 'estimated_richness')
