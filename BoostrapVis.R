#!/usr/bin/env Rscript
# Rscript tax.summary.R file.cons.taxonomy

args = commandArgs(trailingOnly=TRUE)
file = args[1]

taxonomy <- read.csv(file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

tag <- strsplit(file, "[.]")[[1]][1]


if (is.na(args[2])) {
  Rank <- c("Reino", "Filo", "Clase", "Orden", "Familia", "Especie")
} else
  Rank <- c(args[2])

# ==============
## Checking and Load packages ----
# ==============

.cran_packages <- c('dplyr', 'purrr', 'tibble', 'reshape2', 'ggplot2')

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
# Load packages into session, and print package version
sapply(c(.cran_packages), require, character.only = TRUE)

cat("\n.... Ploting!\n")


# x <- read.csv("run012_astr_ASVs.50_COI.wang.taxonomy", header=FALSE, sep="\t", stringsAsFactors=FALSE)
# Ranks <- c("Filo", "Clase", "Orden", "Especie") # BOLD Ranks
# file = "run012_astr_ASVs.50_COI.wang.taxonomy"
# x <- read.csv("run012_astr_ASVs.midori.taxonomy", header=FALSE, sep="\t", stringsAsFactors=FALSE)
# Ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#x <- read.csv("run012_astr_ASVs.IctiosConsenso.wang.tax.summary", header=FALSE, sep="\t", stringsAsFactors=FALSE)
# Ranks <- c("Orden", "Familia", "Especie") # SANGER

tax <- strsplit(taxonomy[,ncol(taxonomy)], ";")
#tax <- sapply(tax, "[", c(1:length(Ranks)+1))
tax <- sapply(tax, "[", c(1:length(Rank)))
tax <- as.data.frame(t(tax))
head(tax)


T <- as.data.frame(apply(tax, 2, 
    function(x) gsub("\\(.*$", "",  x, perl=TRUE)),
     stringsAsFactors = F)

#T <- T[,-c(1)] # For remove midori root Rank
colnames(T) <- Rank
rownames(T) <- taxonomy[,1]

barplot(table(T[,2]), horiz = TRUE, las=2)

Ph <- as.data.frame(table(T[,2]))
spp <- as.data.frame(table(T[,length(Rank)]))


# boostrap
B <- as.data.frame(apply(tax, 2, 
    function(x) gsub("[A-z||()]", "",  x, perl=TRUE)),
     stringsAsFactors = F)

B <- apply(B, 2, as.numeric)
B <- data.frame(B)
# B <- B[,-c(1)]
colnames(B) <- Rank   
rownames(B) <- taxonomy[,1]

boots <- reshape2::melt(B) #, id.vars=colnames(B), variable.name = "Process", value.name = "Bootstrap")
colnames(boots) <- c("Nivel", "Bootstrap")


a <- ggplot(boots, aes(x = Bootstrap))
a <- a +  geom_density(aes(fill = Nivel, y = ..count..), alpha = 0.4) + 
                            scale_fill_brewer(palette = "Dark2") +
                            #scale_color_brewer(palette = "Dark2") +
                            facet_wrap(~Nivel) + 
                            theme_minimal()


write.table(T,
            file=paste0(tag, ".taxonomy"), 
            sep="\t", 
            row.names = TRUE, 
            col.names = TRUE
                   )

 write.table(B,
            file=paste0(tag,".PROBS.taxonomy"), 
            sep="\t", 
            row.names = TRUE, 
            col.names = TRUE
                   )


# Dimentionally of ASVs bootstrap

quit(save = 'no')

pca <- prcomp(t(B), scale=TRUE)
pca <- prcomp(B, scale=TRUE)

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
 
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
 

pca.data <- data.frame(Sample=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2], Filo=T$Filo, Bootstrap=B$Especie)

#T$Sample <- rownames(T)
#T %>%
#    group_by(Sample) %>%
#    inner_join(pca.data, by = "Sample") -> pca.data.plot

ggplot(data=pca.data, aes(x=X, y=Y, size=Bootstrap, color=Filo)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2") +
  #geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_minimal() +
  ggtitle("Bootstrap de Nivels Taxonomicos - PCA")


## ====

quit(save='no')

count <- read.csv("../../run012_astr_ASVs_count.table", sep="\t")
ict <- count[,c(54,55)]

sanger <- data.frame(B, Average=rowSums(B)/ncol(B), BestID=T$Especie, ict, stringsAsFactors=FALSE)

ftest <- sanger[sanger$X04.P68.ICT != 0,] # en funcion de la P114 / P68 hacemos el segundo cutoff
# ftest <- ftest[ftest$Average >=50,] # quitamos malas asignaciones


#aggreg <- aggregate(ftest, list(ftestBestID), sum)


predicted <- as.data.frame(table(ftest$BestID))[,2]
obs <- rep(1, length(predicted))


library(caret) 
xtab <- table(predicted, obs)
confusionMatrix(xtab) # tiene que ser matriz de factores ,no numerica


aggregate(ftest$X04.P68.ICT, by=list(ftest$BestID), FUN=sum)


# ==============
# Outputs config
# ==============

label <- strsplit(list.files(data_path)[1], "-")[[1]]
date <- format(Sys.time(), "%Y%m%d")
out_prefix <- paste0("run",label[1], n_test) # Output prefix (for plots and files)
run <- paste0("run", label[1], "_", date,"_COI")

# ================
# Outputs in `pwd`:
# ================
dir = getwd()
out_path <- file.path(dir, "dada2_asv", run)
system(command = paste0("mkdir -p ", out_path), intern = F)
