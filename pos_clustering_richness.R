# pos-clustering w/ Lulu ----
library(tidyverse)

rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

# devtools::install_github("tobiasgf/lulu") 


load(paste0(dir, "objects.RData"))
source("~/Documents/GitHub/Estadistica_UABC/anova_and_gaussianity.R")

mtd %>% mutate(Tissue = recode_factor(Tissue, Intestine = "Hindgut", Hepatopancreas = "Midgut", Stomach = "Foregut")) -> mtd

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

load(paste0(dir, "lulu_curated_result.Rdata"))

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

# test features diversity ----

source("~/Documents/GitHub/metagenomics/estimate_richness.R")

obj %>%
  group_by(Family) %>% # instead of family, test by asv or other level
  summarise_at(vars(colNames), sum) %>% # or at family level
  select_at(vars(c(colNames))) %>%
  estimate_richness(., measures = c("Observed", "Shannon", "Chao1")) %>%
  as_tibble(rownames = "Index") %>%
  inner_join(mtd) -> diversityDat

library(ggpubr)

diversityDat %>%
  select(-se.chao1, - Index) %>%
  pivot_longer(cols = c("Observed", "Shannon", "Chao1"), names_to = 'alpha') -> eda_df

eda_df %>%
  mutate(Tissue = paste0(Tissue, "-",alpha)) %>%
  rename("x" = value, "g" = Tissue) %>%
  is_parametric()

qqfun <- function(x) {
  x <- x[x > 0]
  qq <- qqnorm(x, plot.it = F) 
  qq %>% as_tibble()
}

zfun <- function(x) {
  x <- x[x > 0]
  z <- c((x - mean(x)) / sd(x))
  # z %>% as_tibble()
  return(z)
}

eda_df %>%
  group_by(alpha, Tissue) %>% # summarise(mean(value))
  mutate(x = value) %>%
  summarise(qqfun(x)) %>%
  mutate(z = zfun(y)) %>%
  mutate(outlier = ifelse(abs(z)>3, TRUE, FALSE)) %>%
  ggplot(aes(x, y, color = Tissue)) +
  geom_point() + 
  ggpubr::stat_cor(aes(group = Tissue), method = "pearson", cor.coef.name = "R", p.accuracy = 0.001) +
  theme_bw(base_family = "GillSans", base_size = 16) +
  labs(x = "Theorical", y = "Diversity Index",
       caption = "Correlation Coefficient (Pearson) between theorical normal distribution (expected) against (observed) distribution in the diversity Indexes\nThe observed distribution fits a linear relatioship to a normal distribution (alfa 5%)") +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
              se = TRUE, na.rm = TRUE, aes(fill = Tissue)) +
  theme(legend.position = "top") +
  facet_wrap(~alpha, scales = 'free_y') -> psave

ggsave(psave, filename = "qqplot_by_asv.png", path = dir, 
       width = 10.5, height = 6.5)

# eda_df  %>% ggqqplot(., x = "value", color = 'Tissue', add.params = list(linetype = "dashed"),
#            conf.int = TRUE) + facet_wrap(~alpha, scales = 'free_y')


library(rstatix)

eda_df %>%
  group_by(alpha) %>%
  rstatix::anova_test(value ~ Tissue) 

# only shannon present difference within the variances of the fractions/tissues

eda_df %>%
  filter(alpha %in% 'Shannon') %>%
  rstatix::anova_test(value ~ Tissue) -> res.aov

# A significant one-way ANOVA is generally followed up by Tukey post-hoc tests to perform multiple pairwise comparisons between groups

eda_df %>%
  filter(alpha %in% 'Shannon') %>%
  group_by(alpha) %>% tukey_hsd(value ~ Tissue) %>%
  add_significance() -> pwc

# if not parametric (the case of asv level)
# a priori
# eda_df %>%
  # group_by(alpha) %>% kruskal_test(value ~ Tissue)
# and posteriori
# eda_df %>%
#   group_by(alpha) %>% wilcox_test(value ~ Tissue) 
# stop here, there's no significant differences at asv level
# 
# The output contains the following columns:
#   
# estimate: estimate of the difference between means of the two groups
# conf.low, conf.high: the lower and the upper end point of the confidence interval at 95% (default)
# p.adj: p-value after adjustment for the multiple comparisons.

pwc <- pwc %>% add_xy_position(x = "Tissue", scales = 'free_y')

eda_df %>%
  # filter(alpha %in% 'Shannon') %>%
  ggboxplot(., x = "Tissue", y = "value",
          bxp.errorbar = T, bxp.errorbar.width = 0.12, outlier.shape = NA) +
  stat_summary(fun = mean, geom="point", shape=20, 
               size = 3, color="red", fill="red") +
  # facet_wrap(~alpha, scales = 'free_y') +
  stat_pvalue_manual(pwc, hide.ns = T) + 
  theme_bw(base_family = "GillSans", base_size = 14)+
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc),
    x = '', y = 'Diversity Index') + 
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1, vjust = 1, size = 10)) -> psave

ggsave(psave, filename = "alpha_Indexes_by_Order.png", path = dir, 
       width = 5.5, height = 4.2)
  


# abort the next bxplot

diversityDat %>%
  pivot_longer(cols = contains(c("Observed", "Shannon", "Chao1"))) %>%
  filter(!name %in% c("se.chao1")) %>%
  ggplot(aes(x = Tissue, y = value)) +
  stat_boxplot(geom ='errorbar', width = 0.15,
               position = position_dodge(0.6)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.6), outlier.shape=NA) +
  stat_summary(fun = mean, geom="point", shape=20, 
               size = 3, color="red", fill="red") +
  coord_cartesian(ylim=c(0,NA)) +
  theme_bw(base_family = "GillSans", base_size = 14) +
  facet_wrap(~name, scales = 'free') +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1, vjust = 1, size = 10)) +
  labs(x = '', y = 'Diversity Index', caption = 'Red dots = mean') -> psave
  
# ggsave(psave, filename = "alpha_Indexes_by_family.png", path = dir, 
#        width = 5.5, height = 4.2)

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
  pivot_wider(names_from = richness, values_from = value) %>%
  ggpaired(x = "group", y = "Observed", 
           order = c("Raw", "Curated"))
  # facet_grid(~ richness, scales = 'free_y')


df_alfa %>% 
  ggplot(aes(x = Tissue, y = value)) +
  facet_grid(richness ~., scales = "free") +
  geom_bar(aes(fill = group), position="dodge", stat="identity") +
  # geom_boxplot(aes(color = group)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1, vjust = 1, size = 10)) +
  ylab("Richness") + xlab("") +
  scale_fill_manual("", values = c("#3182bd", "#de2d26")) -> psave

ggsave(psave, filename = "richness_raw_and_lulu.png", path = dir, 
       width = 5, height = 5)


# stop here!!!!
# are statistically significant different? compare means!!

source("~/Documents/GitHub/Estadistica_UABC/anova_and_gaussianity.R")

df_alfa %>%
  filter(richness %in% "Shannon") %>%
  # group_by(Tissue, Time, group) %>%
  rename("x" = value, "g" = Tissue) %>%
  # mutate(g = paste0(Tissue, group, collapse = "_")) %>%
  is_parametric()

# por tanto, prueba wilcoxcon para dos muestras independientes (tejido)
# The Wilcoxon rank sum test is a non-parametric alternative to the independent two samples t-test for comparing two independent groups of samples, in the situation where the data are not normally distributed

library(rstatix)

df_alfa %>%
  filter(group %in% 'Raw') -> df_alfa

df_alfa %>%
  filter(richness %in% 'Observed') %>%
  # pivot_wider(names_from = richness, values_from = value) %>%
  group_by(richness) %>%
  rstatix::wilcox_test(value ~ Tissue, paired = F, conf.level = 0.95) %>%
  add_significance() -> stat.test

stat.test <- stat.test %>% add_xy_position(x = "Tissue")

df_alfa %>% 
  filter(richness %in% 'Observed') %>%
  ggplot(aes(x = Tissue, y = value)) +
  # facet_grid(richness ~., scales = "free") +
  geom_bar(fill = 'black', position=position_dodge(0.8), stat="identity") +
  # geom_boxplot(aes(color = group)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1, vjust = 1, size = 14)) +
  ylab("Observed") + xlab("") -> p

p + stat_pvalue_manual(stat.test, label = "p.adj.signif", 
                           remove.bracket = F) +
  labs(caption = "A Mann-Whitney U test was evaluated (confidence = 0.95)")-> p

ggsave(p, filename = "Observed_richness.png", path = dir, 
       width = 5, height = 5)

#

df_alfa %>%
  filter(richness %in% 'Shannon') %>%
  # pivot_wider(names_from = richness, values_from = value) %>%
  group_by(richness) %>%
  rstatix::wilcox_test(value ~ Tissue, paired = F, conf.level = 0.95) %>%
  add_significance() -> stat.test

stat.test <- stat.test %>% add_xy_position(x = "Tissue")

df_alfa %>% 
  filter(richness %in% 'Shannon') %>%
  ggplot(aes(x = Tissue, y = value)) +
  # facet_grid(richness ~., scales = "free") +
  geom_bar(fill = 'black', position=position_dodge(0.8), stat="identity") +
  # geom_boxplot(aes(color = group)) +
  theme_bw(base_size = 16, base_family = "GillSans") +
  theme(axis.text.x = element_text(
    angle = 45, hjust = 1, vjust = 1, size = 14)) +
  ylab("Shannon") + xlab("") -> p

p + stat_pvalue_manual(stat.test, label = "p.adj.signif", remove.bracket = T) +
  labs(caption = "A Mann-Whitney U test was evaluated (confidence = 0.95)")-> p

ggsave(p, filename = "Shannon_richness.png", path = dir, 
       width = 5, height = 5)
# my_comparisons <- list( c("Intestine", "Stomach") )
# 
# library(ggpubr)
# 
# psave + 
#   stat_compare_means(
#     # group.by = "group",
#     paired = TRUE, label = "p.signif", label.y.npc = "top")



# R2 and p values of alfa diversity

df_alfa %>%
  filter(group %in% 'Raw') -> df_alfa

df_alfa %>%
  pivot_wider(names_from = richness, values_from = value) %>%
  ggplot(aes(Observed, Shannon)) +
  facet_grid(~Tissue) +
  geom_point() +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", 
                   p.accuracy = 0.001, label.y = 6) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5, 
              se = TRUE, na.rm = TRUE) +
  theme_bw(base_size = 16, base_family = "GillSans") -> psave


ggsave(psave, filename = "richness_raw_cor.png", path = dir, 
       width = 6.4, height = 4.5)




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
