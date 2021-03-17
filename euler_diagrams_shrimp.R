# test wherer or not prevalence groups are equal by tissues
obj %>%  filter(Phylum %in% keepPhyla) %>% select_if(is.double) %>% colSums() %>% sum()
obj%>% select_if(is.double) %>% colSums() %>% sum()
obj %>% 
  filter(Phylum %in% keepPhyla) %>%
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ab') %>%
  filter(ab > 0) %>%
  # select_at(ranks) %>% 
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = ranks) %>% fill(value) %>%
  pivot_wider(names_from = name) %>%
  select(-id) %>%
  left_join(mtd) %>%
  group_by(Tissue,Family) %>%
  summarise(ab = sum(ab)) %>% arrange(desc(ab)) %>%
  ungroup() %>%
  mutate(ab = ifelse(ab > 0, 1, 0)) %>%
  pivot_wider(names_from = Tissue, values_from = ab, values_fill = 0) %>%
  data.frame() -> euler


length(unique(euler$Family)) == nrow(euler)

rownames(euler) <- euler$Family
euler$Family <- NULL

library(UpSetR)
labels <- c("Hindgut", "Midgut", "Foregut")
n <- 3
grid.col <- ggsci::pal_rickandmorty()(n)
names(grid.col) <- labels


png(filename = paste0(dir, "intersect_families.png"),
    width = 780, height = 780, res = 150)
upset(euler, number.angles = 0, point.size = 3.5, line.size = 2, 
      nintersects = 40, 
      queries = list(
        list(query = intersects, params = list("Hindgut"), color = grid.col[1], active = F),
        list(query = intersects, params = list("Midgut"), color = grid.col[2], active = F),
        list(query = intersects, params = list("Foregut"), color = grid.col[3], active = F)),
      mainbar.y.label = "Taxa Intersections", sets.x.label = "Taxa per Tissue", 
      text.scale = c(1.3, 2, 1, 2, 2, 2),
      order.by = "freq")

dev.off()


sum(rowSums(euler) > 1) # continue with comparison per group using those 76 taxa
sum(rowSums(euler) == 1) #make a heatmap of exclusive taxa, 51 exclusive tax

sum(rowSums(euler) > 1)+sum(rowSums(euler) == 1) == nrow(euler)

taxaNames <- rownames(euler)

taxaNames[rowSums(euler) == 1] -> exclusive_taxa
taxaNames[rowSums(euler) > 1] -> intersected_taxa

save(euler, intersected_taxa, exclusive_taxa, file = paste0(dir, "euler_outputs.RData"))

# ordinations


dat <- obj %>% select_at(colNames) %>% data.frame(row.names = obj$`Feature ID`)
tax <- obj %>% select_at(ranks)  %>% data.frame(row.names = obj$`Feature ID`)


ps = phyloseq(otu_table(dat, taxa_are_rows = TRUE), 
                    tax_table(as(tax, 'matrix')), 
                    sample_data(sam)) %>%
  subset_taxa(Family %in% exclusive_taxa) %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%
  prune_samples(sample_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) sqrt(x/ sum(x))) %>%
  aggregate_taxa( "Family")

otu <- abundances(ps)
metadata <- meta(ps)

library(vegan)
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
rda.result <- vegan::rda(t(otu) ~ factor(metadata$Tissue),
                         na.action = na.fail, scale = TRUE)

plot(rda.result)


out.pcoa <- ordinate(ps, method = "PCoA", distance = "bray", na.rm = T)
evals <- out.pcoa$values$Eigenvalues
plot_ordination(ps, out.pcoa, type = "split",
                label = "Phylum", color = "Tissue", justDF = T) -> ordDF
ordDF %>%
  ggplot(aes(Axis.1, Axis.2)) +
  geom_point(aes(color = Tissue, shape = Tissue)) +
  ggrepel::geom_text_repel(data = subset(ordDF, id.type == "Taxa"),
                           aes(label = Phylum), size = 4, max.overlaps = 20) +
  # ggforce::geom_mark_ellipse(aes(group = Tissue)) + # label = Tissue
  facet_grid(~id.type) +
  ggsci::scale_color_startrek()
  
# conviene hacer el RDA y editarlo con extraGrid -----

# set left-panel (Phylum) ordering 

TimeLev <- c("Farm", 0,20,40,60,80)

obj %>% filter(Family %in% exclusive_taxa) %>% 
  mutate_if(is.double, function(x) sqrt(x / sum(x))) %>%  
  pivot_longer(cols = colNames, 
               names_to = "Index", values_to = 'ab') %>%
  filter(ab > 0) %>%
  inner_join(sam) %>%
  mutate(Time = str_replace_all(Time, c("Day" = ""))) %>%
  mutate(Time = factor(Time, levels = TimeLev)) -> dataHeat

dataHeat %>% group_by(Phylum) %>% summarise(t = sum(ab)) %>% arrange(desc(t)) %>% pull(Phylum) -> PhylumLevel

dataHeat %>%
  mutate(Family = factor(Family, levels = exclusive_taxa),
         Phylum = factor(Phylum, levels = PhylumLevel),
         ab = ifelse(ab == 0, NA, ab)) %>%
  ggplot(aes(y = Family, x = Index, fill = ab)) +
  geom_tile(color = "black") +
  ggsci::scale_fill_material("blue-grey",  
                             name="Relative\nAbundance\n", 
                             na.value = 'white'
                             # limits = c(0,100),
                             # labels = scales::percent_format(scale = 1)
                             ) +
  # ggh4x::facet_nested(~ Tissue + Time, scales = "free", space = "free") +
  ggh4x::facet_nested(Phylum + Class ~ Tissue + Time,
                      scales = "free", space = "free" ) +
  labs(x = NULL, y = NULL) + 
  theme_classic(base_family = "GillSans", base_size = 14) +
  guides(fill = guide_colorbar(barheight = unit(9, "in"), 
                               barwidth = unit(0.5, "in"),
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
