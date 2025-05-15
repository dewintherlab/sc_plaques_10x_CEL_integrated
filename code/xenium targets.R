## Make myeloid marker list for Xenium
dir.create("Xenium", showWarnings = F)

## For genreal myeloid markers, use the macrophage archetypes
# Object
integrated.full.seurat

# Archetypes
unique(integrated.full.seurat$archetype)

# Swap archetypes in for idents
integrated.full.arch.seurat         <- integrated.full.seurat
Idents(integrated.full.arch.seurat) <- integrated.full.arch.seurat$archetype
unique(Idents(integrated.full.arch.seurat))

# Myeloid clusters
myeloid.clusters <- levels(unique(Idents(integrated.full.arch.seurat)))
myeloid.clusters <- myeloid.clusters[c(1,2,3,4,10,13,14,19,22)]


## Get markers
integrated.full.seurat.archetype.markers <- FindAllMarkers(object = integrated.full.arch.seurat, min.pct = 0.1, min.diff.pct = 0.25,only.pos = T)

# Save to disk
saveRDS(integrated.full.seurat.archetype.markers, file = "Seurat_Objects/integrated.full.archetype.markers.Rds")
write.better.table(integrated.full.seurat.archetype.markers, file = "Xenium/markers.txt")

### integrated.full.seurat.archetype.markers <- readRDS(file = "Seurat_Objects/integrated.full.archetype.markers.Rds")

# Save the top9 markers per cluster
sep.markers <- integrated.full.seurat.archetype.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)

# plot the top9 markers per myeloid cluster
for(i in myeloid.clusters){
  # Make violin plots
  customVln(features = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = integrated.full.arch.seurat,
                     name      = paste("xenium/",i, " markers.pdf", sep = ""),
                     cols = c(full_set.colors, archetype.colors, "Foamy" = "purple"),
                     width = 20, height = 30, draw.names = F)
}

VlnPlot(object = integrated.full.arch.seurat, features = c("CD14", "CD68", "SPI1", "CEBPA", "ITGAM", "EGR1"), cols = c(full_set.colors, archetype.colors, "Foamy" = "purple"), ncol = 2) * theme(axis.text.x.bottom = element_blank()) + theme(legend.position = "right")
ggsave("Xenium/mac markers.pdf", width = 20, height = 10)

VlnPlot(object = integrated.full.arch.seurat, features = c("FUT4", "CD33", "IL1RN", "CEACAM8", "CSF3R"), cols = c(full_set.colors, archetype.colors, "Foamy" = "purple"), ncol = 2) * theme(axis.text.x.bottom = element_blank())
ggsave("Xenium/neut markers.pdf", width = 10, height = 10)

VlnPlot(object = integrated.full.arch.seurat, features = c("CD27", "NCAM1", "KLRC1", "FCGR3A"), cols = c(full_set.colors, archetype.colors, "Foamy" = "purple"), ncol = 2) * theme(axis.text.x.bottom = element_blank())
ggsave("Xenium/nk markers.pdf", width = 10, height = 10)
