# pos-clustering w/ Lulu ----

rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

# devtools::install_github("tobiasgf/lulu") 


load(paste0(dir, "objects.RData"))

vsearch.file <- list.files(path = dir, pattern = "match_list.txt", full.names = T)
matchlist <- read.table(vsearch.file , header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

feature_tab <- obj %>% select_if(is.double) %>% 
  data.frame(., row.names = obj$`Feature ID`)

# vegan::rarecurve(t(feature_tab))

# # # # # # # # # # #
# Run lulu algorithm
# # # # # # # # # # # 

library(lulu)

min_r <- 1
min_match <- 98
min_cooccurence <- 0.95 

curated_result <- lulu(feature_tab, matchlist, minimum_ratio_type = "min", 
                       minimum_ratio = min_r, 
                       minimum_match = min_match, 
                       minimum_relative_cooccurence = min_cooccurence)

save(curated_result, file = paste0(dir, "lulu_curated_result.Rdata"))

otu_map <- curated_result$otu_map
parent <- otu_map[otu_map$curated == "parent",]

cat("\n...Number of Processed Features are:", length(otu_map$curated), "\n")
cat("Curated in ~real biological Features (parents):",table(otu_map$curated)[2], "\n")
cat("And", table(otu_map$curated)[1], "daughters (merged) Features\n")

curated_table <- curated_result$curated_table

tax <- obj %>% select_at(ranks) %>%
  data.frame(row.names = obj$`Feature ID`) %>%
  filter(row.names(.) %in% rownames(curated_table))

testing <- identical(rownames(tax), rownames(curated_table))
cat("\n.... Taxonomy and count table has been curated:",testing, "\n")

# test beta div

# betadiv ----

library(vegan)

# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html

beta_mds <- function(df) {
  df %>%
    t() %>%
    vegan::vegdist(., method="bray") %>%
    vegan::metaMDS(.) -> mds
  
  mds_data <- as.data.frame(mds$points)
  mds_data$Index <- rownames(mds_data)
  
  return(mds_data)
  
}
curated_table %>%
  beta_mds(.) -> betaLulu


obj %>% # non_curated matrix
  select_at(vars(c(colNames))) %>%
  beta_mds(.) %>%
  left_join(., mtd) %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  geom_point(data = betaLulu, colour = "grey70", alpha = 0.3, size = 2.5) +
  geom_point(aes(fill = Tissue),shape = 21, colour = "grey33") +
  # ggforce::geom_mark_ellipse(aes(label = Tissue, group = Tissue)) +
  # ggsci::scale_color_rickandmorty() +
  ggsci::scale_fill_rickandmorty() +
  labs(caption = "Multidimensional Scaling (NMDS) of\nPairwise beta diversity indexes (Bray-Curtis index)\n ") +
  theme_light(base_size = 16, base_family = "GillSans") +
  theme(legend.position = "top") -> Pbeta

ggsave(Pbeta, filename = "Bray_curtis_NMDS_lulu.png", path = dir, 
       width = 7, height = 5)

# and alpha for lulu
TimeLev <- c("Farm", 0,20,40,60,80)

diversityDat %>% 
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) -> diversityDat

curated_table %>%
  estimate_richness(., measures = c("Observed", "Shannon", "Chao1")) %>%
  as_tibble(rownames = "Index") %>%
  left_join(diversityDat , by = "Index", suffix = c("_Curated", "_Raw")) %>%
  pivot_longer(cols = contains(c("Observed", "Shannon", "Chao1"))) %>%
  separate(col = name, c("richness", "group"), "_") %>%
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) %>%
  filter(richness %in% c("Shannon", "Observed")) -> df_alfa

df_alfa %>% 
  ggplot(aes(x = Tissue, y = value)) +
  facet_grid(richness ~ Time, scales = "free") +
  geom_bar(aes(fill = group), position="dodge", stat="identity") +
  # geom_boxplot(aes(color = group)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1, vjust = 1, size = 10)) +
  ylab("Richness") +
  scale_fill_manual("", values = c("#3182bd", "#de2d26")) -> psave

ggsave(psave, filename = "richness_raw_and_lulu.png", path = dir, 
       width = 7, height = 5)


# stop here!!!!
# are statistically significant different? compare means!!

my_comparisons <- list( c("Intestine", "Stomach") )

library(ggpubr)

psave + 
  stat_compare_means(group.by = "group",paired = TRUE, label = "p.signif", label.y.npc = "top")

source("~/Documents/GitHub/Estadistica_UABC/anova_and_gaussianity.R")

df_alfa %>%
  filter(richness %in% "Shannon") %>%
  # group_by(Tissue, Time, group) %>%
  rename("x" = value, "g" = Index) %>%
  # mutate(g = paste0(Tissue, group, collapse = "_")) %>%
  is_parametric() %>% view()

# por tanto, prueba t para dos muestras dependientes (raw y curated)
# clean NA to unknown


df_alfa %>%
  filter(richness %in% "Shannon") %>%
  group_by(Tissue, group) %>%
  summarise(
    count = n(),
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
  )

df_alfa %>%
  filter(richness %in% "Shannon") %>%
  wilcox.test(value ~ group, data = ., paired = TRUE)

df_alfa %>%
  filter(richness %in% "Shannon") %>%
  pivot_wider(names_from = group) %>%
  mutate(x = Curated-Raw) %>%
  rename("g" = Tissue) -> Dif

is_parametric(Dif)
n <- nrow(Dif)
calculate_r_critical_coeff(n) 
(mean(Dif$x))/(sd(Dif$x)/sqrt(n))

qt(p= 0.975, df= n-1, lower.tail=F)


for(i in ncol(tax):1) {
  # it will replace 'unknown_unclassified' terms for 'Unknown' term
  tax[,i] <- sub("unknown", "Unknown", tax[,i])
  tax[,i] <- sub("NA.*", "Unknown", tax[,i])
  #   Sub taxa with *nclassified with Undetermined term
  tax[,i] <- sub(".*nclassified.*", 
                 "Undetermined", tax[,i], perl=TRUE)
  #   Instead of Undetermined classified , use NA
  #tax[,i][grep("Undetermined|Unknown|NA",
  tax[,i][grep("Undetermined|NA",
               tax[,i], perl=TRUE)] <- NA
}


otu_richness <- colSums(curated_table>0)
total_otu <- nrow(curated_table)
betadiversity <- total_otu[i]/mean(otu_richness)
