## Data visualization examples

## Load packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("Seurat", quietly = TRUE))
  BiocManager::install("Seurat")
if (!require("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2")

## Load objects
scplaque        <- readRDS(file = "scPlaque 10X data.rds")
full_set.colors <- readRDS(file = "Seurat_Objects/full_set.colors.RDS")

## Make some Plots
DimPlot(scplaque, cols = full_set.colors, pt.size = 1.5) 
VlnPlot(scplaque,features = "IL1B", cols = full_set.colors) + NoLegend()
FeaturePlot(scplaque, features = "IL1B", cols = c("grey", "orangered"), pt.size = 2)

# Zoom in on the macs
mac.pops <- as.character(unique(Idents(scplaque)))
mac.pops <- mac.pops[grep("Macrophages", mac.pops)]
VlnPlot(scplaque,features = "IL1B", idents = mac.pops, cols = full_set.colors) + NoLegend()

# Some things you can tweak
VlnPlot(scplaque,features = "IL1B", idents = mac.pops, cols = full_set.colors, pt.size = 0) + theme(legend.position = "right", axis.text.x.bottom = element_blank())

# Save last plot with custom margins for better visibility
ggsave(file = "IL1B in macs.pdf",  width = 20, height = 10)