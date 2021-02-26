# test normalization methods


rm(list = ls())

# lets EDA 

library(tidyverse)
library(phyloseq)

dir <- '~/Documents/Shrimp_Estefany/'

plot_ords <- function(physeq, guide_legend = "left", captions = "") {
  out.pcoa <- ordinate(physeq, method = "PCoA", distance = "bray")
  # out.pcoa <- ordinate(physeq, method = "PCoA", distance = "unifrac", weighted=TRUE)
  evals <- out.pcoa$values$Eigenvalues
  plot_ordination(logt, out.pcoa, type = "samples",
                  color = "Tissue", shape = "Time") + 
    geom_point(size = 4) +
    coord_fixed(sqrt(evals[2] / evals[1])) +
    ggforce::geom_mark_ellipse(aes(group = Tissue)) + # label = Tissue
    ggsci::scale_color_startrek() +
    labs(caption = captions) +
    theme_light(base_size = 16, base_family = "GillSans") +
    theme(legend.position = guide_legend)
}


ps <- readRDS(paste0(dir, "phyloseq.rds"))
ps2 <- transform_sample_counts(ps, function(x) sqrt(x / sum(x)))
psRA <- transform_sample_counts(ps, function(x){x / sum(x)})

psHelli <- psRA
otu_table(psHelli) <- otu_table(vegan::decostand(otu_table(psHelli), method = "hellinger"), taxa_are_rows=TRUE)

library(patchwork)
p1 <- plot_ords(ps2, captions = "sqrt(RA) ords (PCoA - unifran)", guide_legend = "top") 
p2 <- plot_ords(psRA, captions = 'RA ords (PCoA - unifran)', guide_legend = "none")
p3 <- plot_ords(psHelli, captions = 'RA ords (PCoA - unifran)', guide_legend = "none")

p1+p2 + plot_layout(guides='collect', ) & theme(legend.position='top')-> save1

ggsave(save1, filename = "normalization_methods.png", path = dir, 
       width = 16, height = 16)

# therefore select sqrt(RA)  normalization and save
saveRDS(ps2, paste0(dir, "ps_normalized.rds"))
