# http://www.sthda.com/english/articles/32-r-graphics-essentials/132-plot-grouped-data-box-plot-bar-plot-and-more/

# for 18S - ASV ----
# count <- multirun_ASVs_18S_count.table.filtered.csv
# tax <- multirun_ASVs_18S_wang.taxonomy.filtered.csv
# metadata <- multirun_metadata.csv
# path <- /Users/cigom/metagenomics/Franccesco/multirun_xixim_18S 

# for 18S - OTUS (lulu or full data?) ----
# count <- lulu.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.shared
# tax <- lulu.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.cons.taxonomy
# path <- 

# run diversity, ktones distributon and abundant tax group, a


# Set colors
color_Crucero <- c('X04'='#8c510a','X05'= '#bf812d', 'X06'='#dfc27d')
color_strategy <- c('OTU'='#253494', 'ASV'='#636363')

Ranks <- c("Kingdom", "Phylum", "Class", "Order", 
           "Family", "Genus", "Species")

# Set paths ----
p1 <- '/Users/cigom/metagenomics/MG_18S/multirun_downstream_OTUs'
p2 <- '/Users/cigom/metagenomics/Franccesco/multirun_xixim_18S'

out_path <- '~/metagenomics/'

url <- 'https://raw.githubusercontent.com/RJEGR/metagenomics/master/readtx.R'
source(url)


reads_n_features <- function(physeq, ktone = 1, ...) {
  
  psample <- subset_samples(physeq, ...)
  psample <- prune_taxa(taxa_sums(psample) > ktone, psample)
  
  out <- psmelt(psample)
  out <- filter(out, Abundance > 0) %>% select(OTU, Abundance, Crucero, Ranks)
  out <- reshape2::melt(out, 
                       measure.vars = Ranks, 
                       variable.name = 'Rank', value.name='name')
  
  out <- filter(out, name != 'Undetermined_R')
  
  nr <- aggregate(out[,'Abundance'], 
                  by=list(out$Rank, out$Crucero), FUN = sum)
  
  nasv <- function(x) {length(unique(x))}
  
  nf <- aggregate(out[,'OTU'], by=list(out$Rank, out$Crucero), FUN = nasv)
  
  ng <- aggregate(out[,'name'], by=list(out$Rank, out$Crucero), FUN = nasv)
  
  total_reads <- sum(sample_sums(psample)) 
  total_features <- nrow(otu_table(psample))
  
  pct_r <- round(nr[,3] / total_reads * 100, 3) # reads percent
  pct_a <- round(nf[,3] / total_features * 100, 3) # features percent
  # sanity check 
  if(identical(nr[,1], nf[,1])) {
    
    n <- data.frame(Rank = nr[,1], Crucero = nr[,2], 
                    nreads = nr[,3], pct_r = pct_r, nasvs = nf[,3],
                    pct_a = pct_a,
                    ngroup = ng[,3])
  } else
    n <- data.frame(Rank = nr[,1], Crucero = nr[,2], 
                    nreads = nr[,3], pct_r = pct_r, nasvs = nf[,3],
                    pct_a = pct_a,
                    ngroup = ng[,3])
  
  return(n)
  
}
# Set counts tables ----
c1 <- paste0(p1 ,'/','lulu.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.shared')

c2 <- paste0(p2, '/','multirun-All-XIXIMI_Geo-18S-PE_table.biom')

# Set taxa tables ----
t1 <- paste0(p1 ,'/','lulu.trim.contigs.good.unique.pick.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.cons.taxonomy')

t2 <- paste0(p2, '/','multirun_XIXIM_18S_PE_rep_seqs.w2pr2_worms_API02.wang.taxonomy')


# Load data ----
# count table (c)
c_otus <- read.table(c1)

library(biomformat)

c_asvs <- read_biom(c2)
c_asvs <- data.frame(as(biom_data(c_asvs), "matrix"))

# taxonomy (t)
read_rdp2 <- function(file, header = T) {
  
  tax <- read_rdp(file, header = header)
  
  tax <- tax[-8]
  
  colnames(tax) <- c(Ranks, "SL")
  
  tax[is.na(tax)] <- 'Undetermined_R'
  
  tax <-data.frame(row.names = rownames(tax),
                   mutate_all(data.frame(tax), funs(str_replace_all(., c(" "="_")))))
  return(tax)
}

dim(t_otus <- read_rdp2(t1, header = T))
dim(t_asvs <- read_rdp2(t2, header = F))

# plot(table(t_asvs$SL+1))
# plot(table(t_otus$SL))
# Test resolution table ----

r1 <- as.integer(t_otus$SL)+1
r2 <- as.integer(t_asvs$SL)+1

ntax_lev <- rbind(
  data.frame(ntaxa = data.frame(table(r1))[,2],
             pct = data.frame(table(r1))[,2] / sum(data.frame(table(r1))[,2])  * 100,
             Rank =  Ranks[as.integer(names(table(r1)))], 
             Strategy = 'OTU'),
  data.frame(ntaxa = data.frame(table(r2))[,2],
             pct = data.frame(table(r2))[,2] / sum(data.frame(table(r2))[,2])  * 100,
             Rank =  Ranks[as.integer(names(table(r2)))], 
             Strategy = 'ASV'))

ntax_lev <- ntax_lev[!is.na(ntax_lev$Rank),]
ntax_lev$Rank <- factor(tax_lev$Rank, levels = Ranks)

levplot <- ggplot(tax_lev, aes(x = Rank, y = pct, group = Strategy, fill=Strategy, color=Strategy)) + 
  geom_point(size=2, alpha=0.6) + geom_line(size=1, alpha=0.6, linetype="dashed") +
  geom_text(aes(label=ntaxa), size=4, vjust=1, color = 'black') +
  # geom_text(aes(label=diff), size=3, vjust=-0.2, color = 'black') +
  scale_color_manual(values = color_strategy) +
  guides(color = FALSE, fill = FALSE) +
  theme_bw() +
  labs(y = 'Total features (%)', x = '')

ggsave(levplot, path = out_path, filename = 'otus_vs_asvs_resolution_18S.png')


# create metadata ----
# (m) of cruise (c) and stations (s)

m1 <- data.frame(
  row.names = names(c_otus),
  Crucero = sapply(strsplit(names(c_otus), "[_]"), `[`, 2),
  Station =  sapply(strsplit(names(c_otus), "[_]"), `[`, 3),
  Strategy = 'OTU')

m2 <- data.frame(
  row.names = names(c_asvs),
  Crucero = sapply(strsplit(names(c_asvs), "[.]"), `[`, 2),
  Station = sapply(strsplit(names(c_asvs), "[.]"), `[`,1),
  Strategy = 'ASV')

# Estimate richness ----
r_otus <- estimate_richness(c_otus)
r_asvs <- estimate_richness(c_asvs) # include singletos bro!

r_otus <- cbind(r_otus, m1)
r_asvs <- cbind(r_asvs, m2)

r <- rbind(r_otus, r_asvs)

# Plots 
# richness  (r) ----


library(ggplot2)

p <- ggplot(r, aes(x = Crucero, y = Shannon, color = Strategy)) 

box <- p +
  geom_boxplot(notch = TRUE) +
  geom_jitter(shape=19, size=2, alpha = 3/5,
              position=position_jitterdodge(dodge.width=0.7, 
                                            jitter.width=0.2)) +
  stat_summary(
    aes(color = Strategy),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4
  ) +
  scale_color_manual(name=NULL,
                     values=color_strategy) +
  labs(title="",
       subtitle = "",
       x=NULL,
       y="Alfa diversity") +
  theme_classic()

jitter <- p + geom_jitter(
  aes(shape = Strategy, color = Strategy), 
  position = position_jitter(0.2), alpha = 3/5,
  size = 2
) +
  stat_summary(
    aes(color = Strategy),
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4
  ) +
  scale_color_manual(name=NULL,
                     values=color_strategy) +
  labs(title="",
       subtitle = "",
       x=NULL,
       y="Alfa diversity") +
  theme_classic()

ggsave(jitter, path = out, filename = 'otus_vs_asvs_jitter_18S.png')
ggsave(box, path = out, filename = 'otus_vs_asvs_box_18S.png')

#

library(ggpubr)


paired <- ggpaired(r, x = "Strategy", y = "Shannon",
         color = "Crucero", line.color = "gray", line.size = 0.4,
         palette = "jco", facet.by = "Crucero")+
  stat_compare_means(paired = TRUE)

ggsave(paired, path = out, filename = 'otus_vs_asvs_paired_18S.png')

# Test tax using the plot of coverage ----


library(phyloseq)

#identical(sort(names(c_otus)), sort(rownames(m1)))
#identical(rownames(c_otus), rownames(t_otus))

p1 <- phyloseq(otu_table(c_otus, taxa_are_rows = TRUE), 
                    tax_table(as(t_otus, 'matrix')), 
                    sample_data(m1))

p2<- phyloseq(otu_table(c_asvs, taxa_are_rows = TRUE), 
               tax_table(as(t_asvs, 'matrix')), 
               sample_data(m2))

# For otus (1)

n04_1 <- reads_n_features(p1, ktone = 1, Crucero == 'X04')
n05_1 <- reads_n_features(p1, ktone = 1, Crucero == 'X05')
n06_1 <- reads_n_features(p1, ktone = 1, Crucero == 'X06')

n1 <- rbind(n04_1, n05_1, n06_1)

n1 <- cbind(n1,  Strategy = 'OTU')

# For asvs (2)

n04_2 <- reads_n_features(p2, ktone = 1, Crucero == 'X04')
n05_2 <- reads_n_features(p2, ktone = 1, Crucero == 'X05')
n06_2 <- reads_n_features(p2, ktone = 1, Crucero == 'X06')

n2 <- rbind(n04_2, n05_2, n06_2)

n2 <- cbind(n2,  Strategy = 'ASV')

n <- rbind(n1, n2)

n$Rank <- factor(n$Rank, levels = Ranks)

n <- n[!is.na(n$Rank),]


scatter <- ggplot(n, aes(x = Rank, y = ngroup, shape = Crucero, color = pct_r, 
                         group = Crucero)) +
  geom_point(size = 3, alpha = 0.7) + 
  geom_path() +
  scale_color_gradient(name = 'Percent of\nSequences', low="blue", high="red") + 
  labs(caption = 'Features determined as k-tones were removed',
       x= 'Taxonomic Resolution',
       y= "Number of groups") + 
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ Strategy, scales = 'free_y')

ggsave(scatter, filename = "otus_vs_asvs_sequence_cov.png", dpi = 200, path = out_path)

# singletones ----
k1 <- data.table(phy_summary(p1), Strategy = 'OTU')
k2 <- data.table(phy_summary(p2), Strategy = 'ASV')

physeq_summary <- rbind(k1, k2)

physeq_summary$taxAcc <- factor(physeq_summary$taxAcc, levels = Ranks)

physeq_summary[, Filter := "Clean"]
physeq_summary[(logSum >= 5), Filter := "Keep"]

table(physeq_summary$Filter)

physeq_summary$Filter <- factor(physeq_summary$Filter, levels = c('Keep', 'Clean'))

# group sequence sizes by ktone
physeq_summary[, seq_group := ifelse((taxSum >= 1 & taxSum <= 10), "A", "B")]
# physeq_summary[(taxSum >= 1 & taxSum <= 10), seq_group := "A"]
# physeq_summary[(taxSum >= 11 & taxSum <= 99), seq_group := "B"]
# physeq_summary[(taxSum >= 100 & taxSum <= 999), seq_group := "C"]
# physeq_summary[(taxSum >= 1000), seq_group := "D"]


# caption = "Features are grouped according to their number of sequences\nA: sum >= 1 & sum <= 10)\nB: sum >= 11 & sum <= 99)\nC: sum >= 100 & sum <= 999)\nD: sum >= 1000")

logplot <- ggplot(physeq_summary, aes(x = nsamp, y = logSum, 
                                      color = Strategy,
                                      shape = seq_group)) + 
  geom_point(size = 2, alpha = 0.7) + 
  scale_shape_manual(values = c(4, 16)) +
  # scale_color_manual(values = c("#E7B800", "#FC4E07")) + 
  scale_color_manual(values = c("#cccaca", "#253494")) +
  guides(color = FALSE, shape = FALSE) +
  theme_classic() +
  theme(legend.position = "top") + 
  facet_grid(Filter ~., scales = 'free_y', space="free") + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  labs(color = '', x = 'Frequency in samples', y = 'log(Number of sequences)')
       

ggsave(logplot, filename = "otus_vs_asvs_kdist.png", dpi = 200 ,path = out_path)

# path all

library(patchwork)
levplot
save <- (logplot + box ) / (levplot) + scatter
ggsave(save, filename = "otus_vs_asvs_summary.png", dpi = 200 ,path = out_path,
       width = 10, height = 10)
