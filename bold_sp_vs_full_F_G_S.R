# RUN FIRST bold.sp_vs_bold_full.R in order to make the alluv object:

x <- melt(alluv[alluv$x_y !=0, c('ASV', 'sp','full')], id.vars  = c('ASV'),
          variable.name = 'DataBase',
          value.name = 'Rank')

x$Rank <- factor(x$Rank, levels = TL[-1])
x$DataBase <- factor(x$DataBase, levels = levels)

dim(out <- bbold(filter(x, Rank == TL[7:9]), fasta_file = fasta_file, count_tbl = count_tbl)) # 351

# dim(out <- out[out$abund > 1, ]) # 349 # removing singletones

summary(out$abund)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#  2.0    19.0    40.0   162.5    86.0 13551.0 
summary(out$seq_size)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 307.0   313.0   313.0   312.8   313.0   316.0


# Filtar con inner_join los resultados de out, en base a los valores de cambio x_y == 0, 
# ie. quitar del resultado out, valores que esten contenidos en la base x_y == 0 y trabajar entonces con abundancia relativa

dim(equals <- bbold(filter(x_y_rank, x_y == 0), fasta_file = fasta_file, count_tbl = count_tbl)) # 351

summary(equals$abund)
summary(equals$seq_size)

test0 <- melt(out, id.vars = c('ASV', 'seq_size', 'abund'), variable.name = 'DataBase', value.name = 'Rank')
test1 <- melt(equals, id.vars = c('ASV', 'seq_size', 'abund'), variable.name = 'DataBase', value.name = 'Rank')

test0 %>% 
  group_by(Rank) %>%
  anti_join(test1, by = "Rank") %>%
  select(ASV) %>%
  as.data.frame() -> test

table(test$DataBase)

test <- test[!duplicated(test),]

out %>%
  group_by(ASV) %>%
  inner_join(test, by = 'ASV') %>%
  select( -Rank) %>%
  as.data.frame() -> out


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

# 2. 

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

