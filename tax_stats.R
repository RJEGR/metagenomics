#!/usr/bin/env Rscript

# How to run:
# Rscript tax.stats.R file.taxonomy file.fasta

#<!-- El siguiente script toma con archivos de entrada la taxonomía asignada con el módulo rdp de mothur classify.seqs(fasta, taxonomy, cutoff, probs=T), y procesa la taxonomía en dos tablas, (1) los valores de confiabilidad (bootstrap) a lo largo de las asignaciones hechas a las secuencias de amplicones y (2) las asignaciones taxonómicas correspondientes. El script no considera ninguna entidad taxonómica (ej. Bergey’s Taxonomic Outlines) e implementa la etiqueta Rank (ej. Rank_1, Rank_2. Rank_N) para referirse al nivel del las asignaciones. Finalmente, dos visualizaciones presentan la distribución de la confiabilidad en base a los grupos taxonómicos identificados y el tamaño del amplicón. -->

args = commandArgs(trailingOnly=TRUE)

taxonomy.file = args[1]
fasta.file <- args[2] # input fasta to calculate nbases seq and plot vs bootstrap

tag <- strsplit(taxonomy.file, "[.]")[[1]][1]

# ==============
## Checking and Load packages ----
# ==============

.cran_packages <- c('dplyr', 'purrr', 'tibble', 'reshape2', 'ggplot2', 'RColorBrewer')
.bioc_packages <- c("Biostrings", "IRanges")

# 1.
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst], dep=TRUE, repos='http://cran.us.r-project.org')
}
# 2.
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(.bioc_packages[!.inst], ask = F)
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

cat("\n 1. Processing taxonomy... \n")

taxonomy.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

tax.split <- strsplit(taxonomy.obj[, ncol(taxonomy.obj)], ";")
max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) # Using the max rank assignation to names the taxonomy object
taxonomy <- as.data.frame(t(taxonomy))
rownames(taxonomy) <- taxonomy.obj[,1]
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

boots.backup <- boots

# Recorte de algun rank con asignacion root;
boots <- boots[apply(boots, 2, mean) != 100]
boots <- boots[apply(boots, 1, max) <= 100,] #remove rows redundant bootstrap
boots <- boots[apply(boots, 1, min) > 0,] # ""

tax <- tax[names(tax) %in% names(boots)]
tax <- tax[rownames(tax) %in% rownames(boots),]

tax.rank <- 1 # Possition rank in data.frame

rank.stat <- data.frame(lineage = tax[,tax.rank], Boots = boots[,tax.rank])

lineage.stats <- aggregate(rank.stat[,2], by=list(rank.stat[,1]), 
                    FUN = function(x) c(mean = mean(x), 
                                        median = median(x),
                                        sd =sd(x),
                                        ntaxa = length(x),
                                        mintax = min(x),
                                        maxtax = max(x) ))
                                        #quantile(x, probs = c(0.1, 0.25, 0.75, 1)) ) )

lineage.stats <- data.frame(lineage = lineage.stats[,1], lineage.stats[,2])

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

 write.table(lineage.stats,
            file = paste0(tag,".PROBS.taxonomy.stats"), 
            sep=" ", 
            append = FALSE, quote = TRUE,
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

cat('\n Visualizing some data\n')


seqs <- readDNAStringSet(fasta.file)

seqs <- seqs[names(seqs) %in% rownames(boots),]

seqs.width <- width(seqs)


plot <- data.frame(nbases = seqs.width, boots)

plot.melt <- melt(plot, id.vars='nbases', 
                        value.name = 'boots',
                        variable.name = "Rank")

plot.melt$Asignación <- ifelse(plot.melt$boots >= 90, 'Confiable', 'Poco confiable')

# 1.
p <- ggplot(plot.melt, aes(x=nbases, y=boots, color=Rank, group = Rank ))
p +
  geom_point(alpha = 0.4, aes(shape = Asignación)) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~Rank) +
  theme_minimal() +
  labs(title = "Confiabilidad de las asignaciones taxonómicas",
       subtitle = "La confiabilidad a lo largo de la profundidad taxonómica es mostrada en los paneles",
       caption = "Se considera el tamaño del amplicón para determinar el efecto dentro de la asignación", 
       x = "Tamaño del amplicón", y = "% Bootstrap") +
  guides(color = FALSE)

# 2.
colourCount = nrow(lineage.stats)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

p <- NULL
p <- ggplot(lineage.stats, aes(x=sd, y=mean, color=lineage, size = ntaxa ))
p + geom_point() +
scale_color_manual(values = c(getPalette(colourCount)), na.value = "grey", guide = guide_legend(ncol=2)) +
  theme_minimal() +
  labs(title = "Grupos taxonómicos asignados",
       subtitle = "Proporción de asignaciones taxonómicas",
       caption = "La media y desviación estándar (sd) son calculadas para cada grupo taxonómico", 
       x = "% Bootstrap (sd)", y = "% Bootstrap (mean)")

# 3. 
boots.stat <- aggregate(plot.melt[,3], by=list(plot.melt[,2]), 
                    FUN = function(x) c(mean = mean(x), 
                                        median = median(x),
                                        sem =sd(x)/ sqrt(length(x)),
                                        sd = sd(x),
                                        sem_min = mean(x) - 2*sd(x)/ sqrt(length(x)),
                                        sem_max = mean(x) + 2*sd(x)/ sqrt(length(x)),
                                        ntaxa = length(x[x >= 99 ]),
                                        mintax = min(x),
                                        maxtax = max(x) ))

boots.stat <- data.frame(Rank = boots.stat[,1], boots.stat[,2])

pd <- position_dodge(0.1) 

ggplot(boots.stat, aes(x = Rank, y = mean, group = Rank)) + 
    geom_errorbar(aes(ymin=sem_min, ymax=sem_max), width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, aes(size = ntaxa)) + 
    scale_color_brewer(palette = "Blues", direction = -1) + theme_minimal() +
    labs(title = "Confiabilidad de las asignaciones taxonómicas",
       subtitle = "Se muestra la media del valor del bootstrap a lo largo de los niveles taxonómicos",
       caption = "El tamaño de los puntos presenta el número de taxones asignados con un valor igual o mayor a 99 del bootstrap.", 
       x = "Ranking taxonómico", y = "% Bootstrap (Media)")

cat('\n Done! \n')

quit(save='no')

