# RUN FIRST bold.sp_vs_bold_full.R in order to make the alluv object:
# difine input:
# Load data ----

rm(list=ls()); 

full.taxa.obj <- read_rdp(paste0(path_BOLD,'/',bold_all))
sp.taxa.obj <- read_rdp(paste0(path_BOLD,'/',bold_sp))

TL <- c("root","Domain","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TL2 <- c("D", "K", "P", "C", "O", "F", "G", "S")
colnames(full.taxa.obj) <- c(TL, 'SL')
colnames(sp.taxa.obj) <- c(TL, 'SL')

# ranking labels ----

# .... 1
x_ <- NULL
for (i in 1:nrow(sp.taxa.obj)) {
  rl <- sp.taxa.obj$SL[i] + 1
  x_[[i]] <- names(sp.taxa.obj)[rl]
}

# 2.
y_ <- NULL
for (i in 1:nrow(full.taxa.obj)) {
  rl <- full.taxa.obj$SL[i] + 1
  y_[[i]] <- names(full.taxa.obj)[rl]
  
}

# 1.
x_y_rank <- data.frame(ASV = rownames(sp.taxa.obj), 
                       sp= x_, full= y_, 
                       #SL_x = sp.taxa.obj$SL,
                       #SL_y = full.taxa.obj$SL,
                       x_y = sp.taxa.obj$SL - full.taxa.obj$SL,
                       stringsAsFactors = FALSE)

str(input_dat <- x_y_rank[x_y_rank$x_y !=0, c('ASV', 'sp','full')])

# or use alluv instead of input data, alluv[, c('ASV', 'sp','full')]
x <- melt(input_dat, id.vars  = c('ASV'),
          variable.name = 'DataBase',
          value.name = 'Rank')

# x$Rank <- factor(x$Rank, levels = TL[-1])
levels <- c('full', 'sp')

x$DataBase <- factor(x$DataBase, levels = levels)
dim(out <- bbold(filter(x, Rank == TL[7:9]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 363 in filter(x, Rank == TL[2:6]
sum(out$abund) # percent of reads with different assignation , 0.8784525 (less than 1 % of the population)
# dim(out <- out[out$abund > 1, ]) # removing singletones if no relative abundance.

ggplot(out, aes(seq_size, log(abund), color = abund)) +
  geom_point()

# TL[-c(1,7:9)

dim(Phyl_o <- bbold(filter(x, Rank == TL[4]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 464 asvs than reach this level in either, sp or full
dim(Class_o <- bbold(filter(x, Rank == TL[5]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 190 asvs than reach this level in either, sp or full
dim(Ord_o <- bbold(filter(x, Rank == TL[6]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 915 " "
dim(Fam_o <- bbold(filter(x, Rank == TL[7]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 109 " "
dim(Gen_o <- bbold(filter(x, Rank == TL[8]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 226 " "
dim(Sp_o <- bbold(filter(x, Rank == TL[9]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank, rel_ab = TRUE)) # 738 " "
# dim(rest <- bbold(filter(x, Rank == TL[1:3]), fasta_file = fasta_file, count_tbl = count_tbl, x_y_rank = x_y_rank,  rel_ab = TRUE)) # 464 asvs than reach this level in either, sp or full
# bbold_(Fam_o, x_y_rank)

dim(out <- rbind(data.frame(Phyl_o, Rank = TL[4]),
             data.frame(Class_o, Rank = TL[5]),
             data.frame(Ord_o, Rank = TL[6]),
             data.frame(Fam_o, Rank = TL[7]),
             data.frame(Gen_o, Rank = TL[8]),
             data.frame(Sp_o, Rank = TL[9])
))

out$Rank <- factor(out$Rank, levels = TL[4:9])           
sum(out$abund)
ggplot(out, aes(seq_size, log(abund), color = Ref, shape = Ref)) +
  geom_point(alpha = 0.7, aes(size = abund)) + theme_bw() + scale_color_brewer(palette = 'Paired') +
  facet_wrap( ~ Rank, scales = 'free') +
  labs(subtitle = paste0("A set of ", nrow(out), " different ASVs assigned to R rank level due to sequence resolution [x_y != 0] \n",
                         "Normalized Abundance from the total sequence abundance vs sequence size is plotted \n",
                         "The percent of ASVs belong to this set is: ",round(sum(out$abund)), "%"))



# Filtar con inner_join los resultados de out, en base a los valores de cambio x_y == 0, 
# ie. quitar del resultado out, valores que esten contenidos en la base x_y == 0 y 
# trabajar entonces con abundancia relativa

dim(equals <- bbold(filter(x_y_rank, x_y == 0), fasta_file = fasta_file, x_y_rank = x_y_rank, count_tbl = count_tbl)) # 18643
sum(equals$abund) # percent of reads with equal assingation 84.23722

ggplot(equals, aes(seq_size, log(abund))) +
  geom_point(alpha = 0.7, aes(size = abund)) + theme_bw() + scale_color_brewer(palette = 'Paired') +
  labs(subtitle = paste0("A set of ", nrow(out), " different ASVs assigned to R rank level due to sequence resolution [x_y == 0] \n",
                         "Normalized Abundance from the total sequence abundance vs sequence size is plotted \n",
                         "The percent of ASVs belong to this set is: ",round(sum(equals$abund)), "%"))

test0 <- melt(out, id.vars = c('ASV', 'seq_size', 'abund'), variable.name = 'DataBase', value.name = 'Rank')
test1 <- melt(equals, id.vars = c('ASV', 'seq_size', 'abund'), variable.name = 'DataBase', value.name = 'Rank')

test0 %>% 
  group_by(Rank) %>%
  anti_join(test1, by = "Rank") %>%
  select(ASV, Rank) %>%
  as.data.frame() -> test

dim(test) # 76

dim(test <- test[!duplicated(test$ASV),]) # 47

# sanity check
# should print data.frame of 0 columns
ncol(test1[test$ASV %in% test1$ASV]) == 0

if(ncol(test1[test$ASV %in% test1$ASV]) == 0) {
  out %>%
  group_by(ASV) %>%
  inner_join(test, by = 'ASV') %>%
  select( -Rank) %>%
  as.data.frame() -> out
} else {out <- out}

dim(out) == dim(out <- out[!duplicated(out),])


# 0. Select Family
data0 <- rbind(
  data.frame(table(out$sp_F), db = 'sp'),
  data.frame(table(out$full_F), db = 'full')
  
)

names(data0) <- c('lineage', 'Size', 'DataBase')
data0$Rank <- TL[7]

# 1. Select Genus
data1 <- rbind(
  data.frame(table(out$sp_G), db = 'sp'),
  data.frame(table(out$full_G), db = 'full')
  
)

names(data1) <- c('lineage', 'Size', 'DataBase')
data1$Rank <- TL[8]

# 2, Select Order

data2 <- rbind(
  data.frame(table(out$sp_S), db = 'sp'),
  data.frame(table(out$full_S), db = 'full')
  
)

names(data2) <- c('lineage', 'Size', 'DataBase')
data2$Rank <- TL[9]

# 3. Parse results
dim(ntaxa_data <- rbind(data0,data1, data2))
ntaxa_data$Rank <- factor(ntaxa_data$Rank, levels = TL[7:9])
ntaxa_data <- ntaxa_data[order(-ntaxa_data$Size),]

# aggregate colsums by 


tax_sp <- select(out, abund, paste("sp", TL2, sep="_"))
tax_full <- select(out, abund, paste("full", TL2, sep="_"))

# Calculate abundance ----
# coherence with data from ntaxa size
# by


Family <- rbind( data.frame(aglom_ab(tax_sp, 'sp_F'),
                            DataBase = 'sp',
                            Rank = 'Family'),
                 data.frame(aglom_ab(tax_full, 'full_F'),
                            DataBase = 'full',
                            Rank = 'Family'))
# By 
Genus <- rbind( data.frame(aglom_ab(tax_sp, 'sp_G'),
                           DataBase = 'sp',
                           Rank = 'Genus'),
                data.frame(aglom_ab(tax_full, 'full_G'),
                           DataBase = 'full',
                           Rank = 'Genus'))
# By 

Species <- rbind( data.frame(aglom_ab(tax_sp, 'sp_S'),
                           DataBase = 'sp',
                           Rank = 'Species'),
                data.frame(aglom_ab(tax_full, 'full_S'),
                           DataBase = 'full',
                           Rank = 'Species'))
# Plot aglomerated abundance ----

abund_data <- rbind(Family, Genus, Species)
abund_data$Rank <- factor(abund_data$Rank, levels = TL[7:9])
abund_data <- abund_data[order(-abund_data$Size),]


# merge ntaxa and abund
ntaxa_data$abund <- 'nASVs'
abund_data$abund <- "nreads"

data <- rbind(ntaxa_data, abund_data)

data$abund <- factor(data$abund, levels = c('nreads', 'nASVs'))

ggplot(filter(data, abund == 'nASVs'),  aes(lineage, Size, fill = DataBase, group = DataBase)) +
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") + coord_flip() +
  geom_text(aes(label = Size, color = DataBase), size = 3, position = position_dodge(width = 0.5)) +
  scale_color_brewer(palette = "Set1") +
  #facet_grid(Rank~ abund, scales = 'free') +
  facet_wrap(~Rank, scales = 'free') +
  theme_classic(base_size=7) +
  labs(y = 'n ASVs',
    title = 'Number of taxa during database assignation',
    caption = paste0('Using the subset: filter(x, Rank == "Domain","Kingdom","Phylum","Class", "Order") \n Singletones removed'))

# 2. only species 

ggplot(filter(data, Rank == 'Species'),  aes(lineage, Size, fill = DataBase, group = DataBase)) +
  geom_col(width = 0.4, position = position_dodge(), alpha = 0.7) +
  scale_fill_brewer(palette = "Set1") + coord_flip() +
  # geom_text(aes(label = Size, color = DataBase), size = 3, position = position_dodge(width = 0.5)) +
  geom_text(data = subset(data, Rank == 'Species' && abund == 'nASVs'), aes(label=Size,  color = DataBase), 
                position = position_dodge(width = 0.5), check_overlap =FALSE, size = 3) +
  scale_color_brewer(palette = "Set1") +
  facet_grid(Rank~ abund, scales = 'free') +
  theme_bw(base_size=7) +
  labs(y = 'n Size',
       title = 'Number of taxa during database assignation',
       caption = paste0('Using the subset: filter(x, Rank == "Domain","Kingdom","Phylum","Class", "Order") \n Singletones removed'))


# search some rare species assignations
lineage <- c(select(filter(data, abund == 'ntaxa' & DataBase == 'full' & Rank == 'Species'), lineage))

lineage <- lineage$lineage #[2]

out[out$full_S %in% lineage, paste0('sp_', TL2) ]
# df1[with(df1, grepl("B|F", paste(Col2, Col3))),]


# calcular radio de abundancia de ntaxa:nread ----
# 161:17438/161

count.tbl0 <- read.table(paste0(path_BOLD, "/",count_tbl))
# se tienen que contar incluso singletones, sum(count.tbl0[rowSums(count.tbl0) > 1,])

# ternary plot
# see figure 3 from https://www.biorxiv.org/content/biorxiv/early/2019/03/12/575928.full.pdf to replicate it
# install.packages('ggtern')

# install.packages('Ternary')
#library('Ternary')
#TernaryPlot()

library(ggtern)

tern_obj <- subset(x_y_rank_m, x_y != 0  & Rank != 'root')

tern_obj0 <- aggregate(tern_obj[,'DataBase'], by=list(tern_obj[,'Rank']), FUN = table)

tern_shared <- data.frame(table(c(select(subset(x_y_rank_m, x_y == 0 & Rank != 'root'), Rank))))

datavis <- data.frame(Rank = tern_obj0[,1], tern_obj0[,2], shared = tern_shared$Freq[-1])
datavis <- data.frame(Rank = tern_obj0[,1], apply(datavis[-1], 2, function(x) { x / sum(x) * 100})) 

ggtern(data=datavis,aes(full,sp,shared)) + 
  #geom_point() + 
  theme_rgbw() +
  geom_text(aes(label = Rank))

# > datavis
# Rank full  sp shared
# 1  Domain  229 921  11282
# 2 Kingdom  862 616   8890
# 3  Phylum  232 232   1046
# 4   Class  151  39    392
# 5   Order  864  51    898
# 6  Family   72  37    338
# 7   Genus  173  53    632
# 8 Species   52 686  13808

ibrary(ggraph)
library(igraph)
library(tidyverse)
library(viridis)

melt(aggregate(x_y_rank_m[,'DataBase'], by=list(x_y_rank_m[,'Rank']), FUN = table))
data.frame(table(c(select(subset(x_y_rank_m, x_y == 0 & Rank != 'root'), Rank))))

ggraph(x_y_rank_m, layout = 'circlepack', weight="size") + 
  geom_node_circle(aes(fill = as.factor(x_y), color = as.factor(DataBase) )) +
  scale_fill_manual(values=c("0" = "white", "1" = viridis(4)[1], "2" = viridis(4)[2], "3" = viridis(4)[3], "4"=viridis(4)[4])) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black", "3" = "black", "4"="black") ) +
  theme_void() + 
  theme(legend.position="FALSE") 


# okay, tenemos dos posibles hipotesis:
# - a medida que la base de datos aumenta, el algorimo desconoce que hacer con tantas secuencias 
# - a medida que la base de datos aumenta, el algorimo enriquece la informacion para hacer una asignacion mas confiable
# - a medida que la base de datos disminuye, el algoritmo introduce falsos positivos pues asigna la secuencia mas parecida,
# - a medida que la base de datos disminuye, el algoritmo introduce falsos negativos pues no logra asignar el taxon 'unclassify'

# trabajar con el set de secuencias de los asvs unicos (~ 50) y comparar su asignacion entre las dos bases de referencia:
# - evaluar distancia hamming
# - un blast entre las bases de datos y la lista de asvs
# - 
# remover asvs en base a tamano de secuencia,
# codon de paro (a.a), 
# distancia

# seleccionamos asvs ----
# track the list of asvs from databases (full, sp), and use sequences to analyse hamming distance or other
# this task (or other) is neccesary to report false positive or false negative (also define this analysis thorug sanger mock68 )

asv_subset <- out$ASV

seqs0 <- readDNAStringSet(paste0(path_BOLD,'/',fasta_file))
seqs <- seqs0[names(seqs0) %in% asv_subset]
seqs <- seqs[match(asv_subset, names(seqs)),]

# to rename, use last lineage assigned

x_ <- llineage(sp.taxa.obj)
y_ <- llineage(full.taxa.obj)

x_y_lineage <- data.frame(ASV = rownames(sp.taxa.obj), sp= x_, full= y_)
x_y_l_subset <- x_y_lineage[x_y_lineage$ASV %in% asv_subset, ]
x_y_l_subset <- x_y_l_subset[match(asv_subset, x_y_l_subset$ASV),]


library(ShortRead)

# seqs_headers <- paste0(">filtFs_Read", seq(length(seqs)))

#file_name <- "blastn_nr/blastnr_mergers.fasta"
#save_fasta <- c(rbind(seqs_headers, dfM$sequence))
#write(save_fasta, file=file_name)

# dfM$id <- seqs_headers

# for hamming
write(out$ASV, file=paste0(path_BOLD,"/non_shared_asvs.out" ))



taxonomy.file <- paste0(path_BOLD,'/',bold_all)
taxonomy.obj <- read.csv(taxonomy.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
max.rank <- max(lengths(tax.split)) # La funcion lengths esta en R v.3.5
taxonomy <- sapply(tax.split, "[", c(1:max.rank)) 
taxonomy <- as.data.frame(t(taxonomy))
tax0 <- as.data.frame(apply(taxonomy, 2, function(x) gsub("\\(.*$", "",  x, perl=TRUE)), stringsAsFactors = F)

rownames(tax0) <- taxonomy.obj$V1

tax <- unite(tax0, col = 'lineage', sep = ';')

tax_subset <- data.frame(ASV = rownames(tax)[which(rownames(tax) %in% asv_subset)], 
                         lineage = tax[which(rownames(tax) %in% asv_subset), ])

write(tax_subset$lineage, file=paste0(path_BOLD,"/non_shared_taxa.out" ))

#  get the taxonomy, then de db id and sequence from db
# awk 'NR==FNR{a["ASV"$0];next}/^ASV/{f=0;}($0 in a)||f{print;f=1}' non_shared_asvs.out run014_t2_ASVs.BOLD_public_species.wang.taxonomy > test
#  awk '{print $1}' non_shared_asvs.out | xargs -I {} grep "{}" run014_t2_ASVs.BOLD_public_species.wang.taxonomy
# awk 'NR==FNR{a[$1]; next} {for (i in a) if (index($0, i)) print $1}' non_shared_asvs.out run014_t2_ASVs.BOLD_public_species.wang.taxonomy
#### BLAST by hand against nt. Defaults. Megablast.
#### Save hit table (text) (using F2PFFZXG015-Alignment-HitTable.csv)
# and load functions from mock_analysis_hamming.R

seqsM <- as.data.frame(seqs)


allM <- matrix(0, ncol=nrow(seqsM), nrow=1)

colnames(allM) <- rownames(seqsM)
rownames(allM) <- 'abund'

allM['abund',] <- out$abund

dfM <- data.frame(t(allM))

dfM$sequence <- seqsM$x

# comprare vs ref ----
# search in taxonomy.*.csv
dir = '/Users/cigom/metagenomics/db/bold/BOLD_public_trim'
fasta.file = 'dada2_asv/run012_20190329_COI/mock_hits_relax_ASVs.fasta'
reference.file = 'ictio_coi_sanger114.fasta'

ref <- readFasta(reference.file)
refSeqs <- as.character(sread(ref))
strain <- as.character(id(ref))

dfM.vs.ref <- outer(dfM$sequence, refSeqs, evalDist, band=-1)
dfM$refdist <- apply(dfM.vs.ref, 1, min)

# continue with

blast_file <- '/Users/cigom/metagenomics/COI/run012/mock_P68/mock_parameter_definition/blastn_nr/blastnr_mergers-Alignment.txt'
blast_file <- '/Users/cigom/metagenomics/COI/species_resolution_per_db/F2PFFZXG015-Alignment.txt'

dfM$hit <- isHit100(dfM, blast_file)
dfM$oo <- isOneOff(dfM, blast_file)


