load(paste0(dir, "euler_outputs.RData"))
load(paste0(dir, "objects.RData"))

obj %>% select_at(ranks) %>% 
  mutate(id = 1:nrow(obj)) %>%
  pivot_longer(cols = ranks) %>% fill(value) %>%
  pivot_wider(names_from = name) %>%
  select(-id) %>%
  data.frame(row.names = obj$`Feature ID`) -> tax

     
euler %>% as_tibble(rownames = "Family") %>%
  left_join(tax %>% select_at(ranks[2:5]) %>% distinct(Family, .keep_all = T)) %>%
  pivot_longer(cols = names(euler), names_to = "Tissue", values_to = "pa") -> df

df %>% 
  filter(pa > 0) %>%
  filter(Family %in% exclusive_taxa) %>%
  with(., table(Phylum, Tissue)) %>% t()
  
df %>% 
  filter(pa > 0) %>%
  filter(Family %in% intersected_taxa) %>%
  with(., table(Phylum, Tissue)) %>% t()

# df %>% filter(Family %in% exclusive_taxa) -> df
# df %>% filter(Family %in% intersected_taxa) -> df

library(circlize)

set.seed(260220)

n <- length(unique(df$Tissue))
grid.col <- ggsci::pal_rickandmorty()(n)
names(grid.col) <- c("Hindgut", "Midgut", "Foregut")# sort(unique(df$Tissue),decreasing = T)

# establecer paleta de colores de filos

df %>% pull(Phylum) %>% sort() %>% unique() -> labels 

colourCount = length(labels)

library(ggsci)

if(colourCount > 7) {
  getPalette <- colorRampPalette(pal_locuszoom()(7))(colourCount)
} else
  getPalette <- pal_locuszoom()(colourCount)  

names(getPalette) <- labels


filename = paste0(dir, "/intersected_taxa_chord.png")
# intercalate the filename
# paste0(dir, "/exclusive_taxa_chord.png")
width= 1500;height=1500;res = 150
png(filename, width = width, height = height, res = res)


circos.clear()
circos.par(start.degree = 0, gap.degree = 4, 
           track.margin = c(-0.01, 0.01), 
           points.overflow.warning = FALSE)
df %>% 
  filter(pa > 0) %>%
  select(Tissue, Phylum) %>%
  with(., table(Phylum, Tissue)) %>%
  chordDiagram(grid.col = c(grid.col, getPalette),
               directional = -1,
               # annotationTrack = "grid", 
               # direction.type = c("arrows", "diffHeight"),
               preAllocateTracks = 1,
               small.gap = 10)

dev.off()
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1],
#               CELL_META$sector.index,
#               facing = "clockwise",
#               niceFacing = TRUE, adj = c(0, 0.5))
# }, bg.border = NA)

# abline(h = 0.05, lty = 2, col = "#00000080")

# later test some of this https://www.nature.com/articles/s41522-018-0077-y#Sec7
