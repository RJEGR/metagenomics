

rm(list = ls())

dir <- '~/Documents/Shrimp_Estefany/'

load(paste0(dir, "objects.RData"))
# Analysis composition of microbiomes w/ ANCOM-BC ----

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ANCOMBC")

# http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
#
library(ANCOMBC)
library(phyloseq)
library(microbiome)


sam <- data.frame(mtd, row.names = mtd$Index) %>% 
  arrange(match(Index, colNames)) %>% mutate_if(is.character, as.factor)

dat <- obj %>% select_at(colNames) %>% data.frame(row.names = obj$`Feature ID`)
tax <- obj %>% select_at(ranks)  %>% data.frame(row.names = obj$`Feature ID`)
identical(names(dat), rownames(sam))
identical(rownames(dat), rownames(tax))

# and parse
phyloseq = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')), 
                    sample_data(sam)) %>%
  subset_taxa(Family %in% tax_hclust) %>%
  prune_taxa(taxa_sums(.) > 0, .)

ancombc_data <- aggregate_taxa(phyloseq, "Family")

# 

# Run ancombc function

out <- ANCOMBC::ancombc(phyloseq = ancombc_data, formula = "Tissue",
                        p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000,
                        group = "Tissue", struc_zero = TRUE, neg_lb = TRUE,
                        tol = 1e-5, max_iter = 100, conserve = TRUE,
                        alpha = 0.05, global = TRUE)

res = out$res


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

table(df_fig3$star)

# set the contrast (look at the contrasts)
mtd %>%
  with(., table(Tissue, Time))
# colnames(df_fig2)[-1] = paste0(colnames(df_fig2)[-1], "SD")

df_fig1 %>% left_join(df_fig2) %>% left_join(df_fig3) %>%
  filter(logFC != 0) %>% arrange(desc(logFC)) %>%
  mutate(group = str_replace_all(group, c("Tissue" = ""))) %>%
  mutate(wrap = group) %>%
  mutate(group = ifelse(logFC > 0, "g1", "g2")) -> df_fig

df_fig$taxon_id = factor(df_fig$taxon_id, levels = unique(df_fig$taxon_id))

# Data are represented by effect size (log fold change) and 95% confidence interval bars (two-sided; Bonferroni adjusted) derived from the ANCOM-BC model.\nAll effect sizes with adjusted p<0.05 are indicated:

p = ggplot(data = df_fig, 
           aes(x = taxon_id, y = logFC, fill = group, color = group)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = logFC - SE, 
                    ymax = logFC + SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  geom_text(aes(y = logFC+2.5*sign(logFC), label=star), 
            vjust=.7, color="black", position=position_dodge(width = .5)) +
  scale_fill_discrete("Intercept") +
  guides(color = "none") +
  labs(x = NULL, y = "Log fold change", 
       title = "Analysis composition of microbiomes w/ \nbias correction (ANCOM-BC)\n",
       caption = "g1 (ie. Tissue in the facet) and g2 (Hepatopancreas as intercept)\n *significant at 5% level of significance;\n**significant at 1% level of significance;\n***significant at 0.1% level of significance") + 
  theme_classic(base_size = 16, base_family = "GillSans") + 
  facet_wrap(~wrap) +
  coord_flip() +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0),
        plot.caption = element_text(hjust = 0),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(p, filename = "ANCOMBC.png", path = dir, 
       width = 8, height = 6)

# re-doing heatmap

diff_ab_list <- as.character(df_fig$taxon_id)
agglom_lev <- "Family"
TimeLev <- c("Farm", 0,20,40,60,80)

ancombc_data %>%
  subset_taxa(Family %in% diff_ab_list) %>%
  psmelt() %>%
  rename("taxon_id" = OTU, "ab" = Abundance) %>% 
  group_by(Index, Tissue, Time) %>%
  mutate(ra = ab / sum(ab) * 100) %>%
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) -> dataHeat

# sanity check ----

dataHeat %>% summarise(sum(ra))

# set left-panel (Phylum) ordering 
dataHeat %>% group_by(Phylum) %>% summarise(t = sum(ra)) %>% arrange(desc(t)) %>% pull(Phylum) -> PhylumLevel

dataHeat %>%
  mutate(taxon_id = factor(taxon_id, levels = levels(df_fig$taxon_id)),
         Phylum = factor(Phylum, levels = PhylumLevel),
         ra = ifelse(ra == 0, NA, ra)) %>%
  ggplot(aes(y = taxon_id, x = Index, fill = ra)) +
  geom_tile() +
  # facet_grid( ~ as.factor(Tissue) , scales = "free", space = "free", switch = "x") +
  geom_tile(color = "black") +
  ggsci::scale_fill_material("blue-grey",  
                             name="Relative\nAbundance\n", 
                             na.value = 'white',
                             limits = c(0,100),
                             labels = scales::percent_format(scale = 1)) +
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
    panel.spacing = unit(0.007, "lines")) -> heatPlot

ggsave(heatPlot, filename = "heatmap_ANCOM_list.png", path = dir, 
       width = 16, height = 8)
