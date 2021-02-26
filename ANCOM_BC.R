

rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

load(paste0(dir, "objects.RData"))
load(paste0(dir, "euler_outputs.RData"))
# Analysis composition of microbiomes w/ ANCOM-BC ----

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ANCOMBC")

# http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#
library(ANCOMBC)
library(phyloseq)
library(microbiome)
library(tidyverse)

sam <- data.frame(mtd, row.names = mtd$Index) %>% 
  arrange(match(Index, colNames)) %>% mutate_if(is.character, as.factor)

# foregut (stomach), midgut (hepatopancreas), and hindgut (intestine), 

sam %>% mutate(Tissue = recode_factor(Tissue, Intestine = "Hindgut", Hepatopancreas = "Midgut", Stomach = "Foregut")) -> sam

TimeLev <- c("Farm", 0,20,40,60,80)

sam %>% mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) -> sam
# obj %>% filter(Phylum %in% keepPhyla) %>% pull(Class) %>% unique() 

# obj %>% filter(Phylum %in% keepPhyla) %>% 
  # mutate(Family = ifelse(is.na(Family), paste0(Class, "_Incomplete"), Family)) -> obj

# after test wherer or not prevalence groups are equal by tissues
# select the intersected groups

obj %>% select_at(ranks) %>% 
  mutate(id = 1:nrow(obj)) %>%
  pivot_longer(cols = ranks) %>% fill(value) %>%
  pivot_wider(names_from = name) %>%
  select(-id) %>%
  data.frame(row.names = obj$`Feature ID`) -> tax

dat <- obj %>% select_at(colNames) %>% data.frame(row.names = obj$`Feature ID`)
# tax <- obj %>% select_at(ranks)  %>% data.frame(row.names = obj$`Feature ID`)

identical(names(dat), rownames(sam))
identical(rownames(dat), rownames(tax))

# and parse
phyloseq = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')), 
                    sample_data(sam)) %>%
  subset_taxa(Family %in% intersected_taxa) %>%
  prune_taxa(taxa_sums(.) > 0, .)
  # transform_sample_counts(function(x) sqrt(x / sum(x)))


ranks

ancombc_data1 <- aggregate_taxa(phyloseq, "Family")

#  min(sample_sums(ancombc_data1))

# Run ancombc function

# set the contrast (look at the contrasts)
sam %>%
  with(., table(Tissue, Time)) %>% rowSums()

formula <- "Tissue"

out1 <- ANCOMBC::ancombc(phyloseq = ancombc_data1, formula = formula,
                        p_adj_method = "holm", zero_cut = 0.98, lib_cut = 100,
                        group =  formula, struc_zero = TRUE, neg_lb = TRUE,
                        tol = 1e-5, max_iter = 100, conserve = TRUE,
                        alpha = 0.05, global = FALSE)

out_df <- function(out_ancombc) {
  
  res = out_ancombc$res
  
  df_fig1 = data.frame(res$beta * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>%
    pivot_longer(-taxon_id, names_to = "group", values_to = "logFC")
  
  df_fig2 = data.frame(res$se * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>% 
    pivot_longer(-taxon_id, names_to = "group", values_to = "SE")
  
  df_fig3 = data.frame(res$q_val * res$diff_abn, check.names = FALSE) %>% 
    rownames_to_column("taxon_id") %>%
    pivot_longer(-taxon_id, names_to = "group", values_to = "q_val") %>%
    mutate(star = ifelse(q_val <.001, "***", 
                         ifelse(q_val <.01, "**",
                                ifelse(q_val <.05, "*", ""))))
  
  # table(df_fig3$star)
  # colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")
  
  df_fig1 %>% left_join(df_fig2) %>% left_join(df_fig3)
}

out_df(out1) %>%
  filter(logFC != 0) %>% arrange(desc(logFC)) %>%
  mutate(group = str_replace_all(group, c("Tissue" = ""))) %>%
  mutate(wrap = paste0("Hindgut-", group)) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig

# df_fig$taxon_id = factor(df_fig$taxon_id, levels = unique(df_fig$taxon_id))
table(df_fig$wrap)

# compare third group

phyloseq %>%
  subset_samples(Tissue != "Hindgut") %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  aggregate_taxa(., "Family") -> ancombc_data2

# formula <- "Time + Tissue"
formula <- "Tissue"

out2 <- ANCOMBC::ancombc(phyloseq = ancombc_data2, formula = formula,
                        p_adj_method = "holm", zero_cut = 0.98, lib_cut = 100,
                        group =  formula, struc_zero = TRUE, neg_lb = TRUE,
                        tol = 1e-5, max_iter = 100, conserve = TRUE,
                        alpha = 0.05, global = FALSE)
out_df(out2) %>%
  filter(logFC != 0) %>% arrange(desc(logFC)) %>% 
  mutate(group = str_replace_all(group, c("Tissue" = ""))) %>%
  mutate(wrap = paste0(group, "-Midgut")) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig_F_M

table(df_fig_F_M$wrap)
table(df_fig$wrap)
# Data are represented by effect size (log fold change) and 95% confidence interval bars (two-sided; Bonferroni adjusted) derived from the ANCOM-BC model.\nAll effect sizes with adjusted p<0.05 are indicated:

rbind(df_fig_F_M, df_fig) %>%
  ggplot(data = ., 
           aes(x = taxon_id, y = logFC, fill = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = logFC - SE, 
                    ymax = logFC + SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = logFC+2.5*sign(logFC), label=star), 
            vjust=.7, color="black", position=position_dodge(width = .5)) +
  ggsci::scale_fill_uchicago() +
  # scale_fill_discrete("Intercept") +
  # ggsci::scale_fill_gsea(name = "P value (adjusted) ") +
  guides(color = "none") +
  labs(x = NULL, y = "Log fold change", 
       title = "Analysis composition of microbiomes w/ \nbias correction (ANCOM-BC)\n",
       caption = "\n *significant at 5% level of significance;\n**significant at 1% level of significance;\n***significant at 0.1% level of significance") + 
  theme_classic(base_size = 16, base_family = "GillSans") + 
  facet_wrap(~wrap) +
  coord_flip() +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1)) -> p

ggsave(p, filename = "ANCOMBC_all_contrast.png", path = dir, 
       width = 8, height = 8)


# Step 2: adjust the log observed abundances by subtracting the estimated  sampling fraction from log observed abundances of each sample.
adj_ab <- function(out_ancombc, phyloseq) {
  samp_frac = out_ancombc$samp_frac
  # Replace NA with 0
  samp_frac[is.na(samp_frac)] = 0 
  
  # Add pesudo-count (1) to avoid taking the log of 0
  log_obs_abn = log(abundances(phyloseq) + 1) 
  # Adjust the log observed abundances
  log_obs_abn_adj = t(t(log_obs_abn) - samp_frac)
  
  log_obs_abn_adj %>% as_tibble(rownames = "Family") %>% pivot_longer(-Family, names_to = "Index")
}

adj_ab(out1,ancombc_data1) %>%
  right_join(adj_ab(out2, ancombc_data2), by = c("Family", "Index"))

# re-doing heatmap
ancombc_res <- rbind(df_fig_F_M, df_fig) %>% arrange(taxon_id)

diff_ab_list <- ancombc_res %>% pull(taxon_id) %>% unique()

# intersected_taxa

phyloseq %>%
  aggregate_taxa(., "Family") %>%
  # transform_sample_counts(function(x) sqrt(x / sum(x))) %>%
  subset_taxa(Family %in% diff_ab_list) %>%
  psmelt() %>%
  rename("taxon_id" = OTU, "ab" = Abundance) %>% 
  as_tibble() -> dataHeat

# sanity check

phyloseq %>%
  subset_taxa(Family == "Rhodobacteraceae") %>%
  prune_taxa(taxa_sums(.) > 0, .)

dataHeat %>%
  filter(taxon_id %in% "Rhodobacteraceae")


dataHeat %>%
  group_by(taxon_id, Tissue) %>%
  summarise(Prev = sum(ab > 0), Tab = sum(ab)) %>%
  # mutate(Tab = log(Tab + 1)) %>%
  left_join(dataHeat %>% select_at(c("taxon_id", ranks[1:4])) %>% distinct(.keep_all = T)) %>%
  arrange(desc(Tab)) %>%
  group_by(Tissue) %>%
  mutate(Rank = rank(Tab)) -> dataViz

# labels <- dataViz %>% pull(Phylum) %>% unique() %>%
dataViz %>% group_by(Phylum) %>% summarise(t = sum(Tab)) %>% arrange(desc(t)) %>% pull(Phylum) -> labels 

colourCount = length(labels)

library(ggsci)

if(colourCount > 7) {
  getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
} else
  getPalette <- pal_locuszoom()(colourCount)  



# dataViz %>%
#   mutate(Phylum = factor(Phylum,  levels = labels)) %>%
#   mutate(taxon_id = forcats::fct_reorder(taxon_id, Rank)) %>%
#   ggplot() +
#   geom_col(aes(x = taxon_id ,
#                y = Tab, fill = Phylum)) +
#   labs(y = "") +
#   scale_fill_manual(values = getPalette) +
#   facet_grid(Phylum ~ Tissue , scales = "free", space = "free_y") +
#   coord_flip() +
#   theme_classic(base_size = 16, base_family = "GillSans") +
#   theme(strip.background.y = element_blank(), strip.text.y = element_blank())

dataViz %>%
  ggplot(aes(Tab, Prev)) +
  geom_point(aes(color = Tissue)) +
  scale_x_log10() +
  geom_smooth(aes(color = Tissue), se = F) 
  # facet_wrap(~Phylum)
  
dataViz %>%
  arrange(taxon_id) %>%
  group_by(taxon_id) %>%
  mutate(ra = Tab / sum(Tab)) %>%
  mutate(Phylum = factor(Phylum,  levels = labels)) %>%
  mutate(taxon_id = forcats::fct_reorder(taxon_id, Rank)) %>%
  ggplot() +
  geom_col(aes(x = taxon_id ,
               y = ra, fill = Tissue)) +
  ggsci::scale_fill_rickandmorty() +
  # facet_grid(Phylum ~ ., scales = "free", space = "free") +
  ggh4x::facet_nested(Phylum  ~., scales = "free", space = "free", switch = "y") +
  labs(x = NULL, y = NULL) +
  coord_flip() +
  theme_classic(base_size = 18, base_family = "GillSans") +
  theme(strip.background.y = element_blank(),
        strip.text = element_text(margin = margin(1, 1, 1, 1)),
        axis.text.y.left = element_text(angle = 0, hjust = 1, vjust = 1, size = 18),
        axis.text.y.right = element_text(angle = 0, hjust = 1, vjust = 1, size = 12)) -> p2


ggsave(p2, filename = "bar_ANCOM_list.png", path = dir, 
       width = 8, height = 8)

# sanity check ----

dataHeat %>% group_by(Tissue) %>%summarise(sum(ab))

# set left-panel (Phylum) ordering 
dataHeat %>% group_by(Phylum) %>% summarise(t = sum(ab)) %>% arrange(desc(t)) %>% pull(Phylum) -> PhylumLevel

dataHeat %>%
  mutate(Phylum = factor(Phylum, levels = PhylumLevel)) %>%
  mutate(ab = log(ab + 1)) %>%
  ggplot(aes(y = taxon_id, x = Index, fill = ab)) +
  geom_tile() +
  # facet_grid( ~ as.factor(Tissue) , scales = "free", space = "free", switch = "x") +
  geom_tile(color = "black") +
  ggsci::scale_fill_material("blue-grey",  
                             name="Relative\nAbundance\n", 
                             na.value = 'white') +
                             # limits = c(0,100),
                             # labels = scales::percent_format(scale = 1)) +
  # ggh4x::facet_nested(~ Tissue + Time, scales = "free", space = "free") +
  ggh4x::facet_nested(~ Tissue + Time,
                      scales = "free", space = "free" ) +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_family = "GillSans", base_size = 14) +
  guides(fill = guide_colorbar(barheight = unit(6, "in"), 
                               barwidth = unit(0.3, "in"),
                               ticks.colour = "black", 
                               frame.colour = "black",
                               label.theme = element_text(size = 14))
  ) +
  
  theme(
    title = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1, size = 12),
    strip.text.y = element_text(
      angle = 0, 
      size = 14),
    strip.background = element_rect(colour = "black", 
                                    fill = "transparent",
                                    size = 0.4),
    panel.spacing = unit(0.007, "lines"))

ggsave(heatPlot, filename = "heatmap_ANCOM_list.png", path = dir, 
       width = 16, height = 8)
