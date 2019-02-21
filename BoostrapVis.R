#!/usr/bin/env Rscript
# Rscript tax.stats.R file.cons.taxonomy

args = commandArgs(trailingOnly=TRUE)
taxonomy.file = args[1]

tag <- strsplit(taxonomy.file, "[.]")[[1]][1]

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

cat("\n 1. Processing taxonomy... \n")

taxonomy.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) # Using the max rank assignation to names the taxonomy object
taxonomy <- as.data.frame(t(taxonomy))
tax <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)


rank.names <- vector(max.rank, mode="character")

for (i in 1:max.rank) {
  rank.names[i] <- paste("Rank", i, sep="_")
}

colnames(tax) <- rank.names
rownames(tax) <- taxonomy.obj[,1]

# cat("\n.... Ploting!\n")
#Ph <- as.data.frame(table(T[,2]))
#spp <- as.data.frame(table(T[,length(Rank)]))


cat("\n 2. Bootstrap stats ... \n")

boots <- as.data.frame(apply(taxonomy, 2, function(x) gsub("[A-z||()]", "",  x, perl=TRUE)), stringsAsFactors = F)
boots <- apply(boots, 2, as.numeric)
boots <- data.frame(boots)

boots[is.na(boots)] <- 0

colnames(boots) <- rank.names   
rownames(boots) <- taxonomy.obj[,1]

# Recorte de algun rank con asignacion root;
boots <- boots[apply(boots, 2, mean) != 100]
tax <- tax[names(tax) %in% names(boots)]

tax.rank <- 2 # Possition rank in data.frame

rank.stat <- data.frame(rank = tax[,tax.rank], Boots = boots[,tax.rank])

# # # Testing
ggplot(rank.stat, aes(Boots, ..count.., colour=rank, fill=rank)) +
    geom_histogram(alpha=0.35, position = "stack", bins = 100) +
    scale_fill_brewer(palette = "Paired" ) +
    scale_color_brewer(palette ="Paired" ) +
    facet_wrap(~ rank, scales="free_y") + theme_minimal()  
# # # 

boots.stat <- aggregate(rank.stat[,2], by=list(rank.stat[,1]), 
                    FUN = function(x) c(mean = mean(x), 
                                        median = median(x),
                                        sd =sd(x),
                                        ntaxa = length(x),
                                        mintax = min(x),
                                        maxtax = max(x) ))
                                        #quantile(x, probs = c(0.1, 0.25, 0.75, 1)) ) )

boots.stat <- data.frame(Rank = boots.stat[,1], boots.stat[,2])

cat("\n 3. Writing results ... \n")

# 1.
write.table(tax,
            file = paste0(tag, ".taxonomy"), 
            sep=" ",
            append = FALSE, quote = FALSE,
            row.names = TRUE, 
            col.names = TRUE
                   )

# 2.

 write.table(boots.stat,
            file = paste0(tag,".PROBS.taxonomy.stats"), 
            sep=" ", 
            append = FALSE, quote = FALSE,
            row.names = FALSE, 
            col.names = TRUE
                   )

# 3.                  

 write.table(boots,
            file = paste0(tag,".PROBS.taxonomy"), 
            sep=" ", 
            append = FALSE, quote = FALSE,
            row.names = TRUE, 
            col.names = TRUE
                   )

cat('\n Done! \n')

quit(save='no')

# input fasta better to further analysis and, 
# calculate nbases seq from results and plot as follow:

nbases <- read.csv('run012_relax_ASVs.summary', header=TRUE, sep="\t", stringsAsFactors=FALSE)

nbases <- nbases['nbases']

id.vars=colnames(B), variable.name = "Process", value.name = "Bootstrap")
colnames(boots) <- c("Nivel", "Bootstrap")


plot <- data.frame(nbases = nbases, boots)
plot.melt <- melt(plot, id.vars='nbases', 
                        value.name = 'Bootstrap',
                        variable.name = "Rank")
# queda bien, pero hace falta describir que esta indicando
p <- ggplot(plot.melt, aes(x=nbases, y=Bootstrap, color=Rank, group = Rank ))
p +
  geom_point(alpha = 0.4) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~Rank) +
  theme_minimal()


a <- ggplot(datavis, aes(ntaxa)) +
    geom_density(aes(fill = Rank, y = ..count..), alpha = 0.4) + 
                            scale_fill_brewer(palette = "Dark2") +
                            #scale_color_brewer(palette = "Dark2") +
                            facet_wrap(~Nivel) + 
                            theme_minimal()

ggplot(datavis, aes(x=Rank, y=ntaxa, fill=Rank)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=ntaxa, ymax=ntaxa+sd), width=.2,
                 position=position_dodge(.9))



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
