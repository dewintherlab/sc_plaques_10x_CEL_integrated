#=================================================================================================================================
## Comparison of our data with wirka et al.
## Specifically, can we map our populations and do we see the same trajectories.
##================================================================================================================================
dir.create("wirka_validation_results", showWarnings = F)

## Load the wirka data
# Load the object into the env
wirka.seurat <- readRDS("raw_data/Wirka_2019_plaqview.rds")

# Clean up
wirka.seurat <- UpdateSeuratObject(wirka.seurat)
wirka.seurat <- RunPCA(wirka.seurat)
wirka.seurat <- RunUMAP(wirka.seurat, dims = 1:30, return.model = T)

# Check it out
wirka.seurat
unique(Idents(wirka.seurat))
row.names(GetAssayData(wirka.seurat))

# Plot the given populations
DimPlot(wirka.seurat)
customUMAP(object = wirka.seurat, pt.size = 1.5, label = F, title = "wirka et al. data", shuffle = T, seed = 666, file.name = "wirka_validation_results/wirka.UMAP.pdf")


## Transfer our labels by reference mapping
# Find anchors
wirka.anchors <- FindTransferAnchors(reference         = integrated.full.seurat,
                                   query               = wirka.seurat,
                                   dims                = 1:30, 
                                   reference.reduction = "pca")

# Predict labels
predictions <- TransferData(anchorset = wirka.anchors, 
                            refdata   = Idents(integrated.full.seurat),
                            dims      = 1:30)

# Add to the wirka seurat object
wirka.seurat <- AddMetaData(wirka.seurat, metadata = predictions)

# Map to ref UMAP
integrated.full.seurat <- RunUMAP(integrated.full.seurat, dims = 1:50, return.model = T)
wirka.seurat <- MapQuery(anchorset           = wirka.anchors, 
                       reference           = integrated.full.seurat, 
                       query               = wirka.seurat,
                       refdata             = Idents(integrated.full.seurat), 
                       reference.reduction = "pca", 
                       reduction.model     = "umap")

# Visualize
customUMAP(object = wirka.seurat, group.by = "predicted.id", pt.size = 1.5, label = F, title = "wirka et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "wirka_validation_results/wirka.full_ref_map.UMAP.pdf", cols = full_set.colors, plot.width = 15)
customUMAP(object = wirka.seurat, group.by = "predicted.id", reduction = "ref.umap", pt.size = 1.5, label = F, title = "wirka et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "wirka_validation_results/wirka.full_ref_umap.UMAP.pdf", legend.pos = "right", cols = full_set.colors, plot.width = 15)

## Now subset the macropahges and try again
# Subset the wirka object
mye.wirka.seurat <- subset(wirka.seurat, idents = "Mø")

# Re-center
mye.wirka.seurat <- RunPCA(mye.wirka.seurat)
mye.wirka.seurat <- RunUMAP(mye.wirka.seurat, dims = 1:30, return.model = T)

# Check it out
mye.wirka.seurat
unique(Idents(mye.wirka.seurat))

# Plot the given populations
DimPlot(mye.wirka.seurat)
customUMAP(object = mye.wirka.seurat, pt.size = 1.5, label = F, title = "wirka et al. data", shuffle = T, seed = 666, file.name = "wirka_validation_results/wirka.mye.UMAP.pdf")


## Transfer our labels by reference mapping
# Find anchors
DefaultAssay(integrated.mye.seurat) <- "integrated"
mye.wirka.anchors <- FindTransferAnchors(reference           = integrated.mye.seurat,
                                       query               = mye.wirka.seurat,
                                       dims                = 1:50, 
                                       reference.reduction = "pca")

# Predict labels
predictions <- TransferData(anchorset = mye.wirka.anchors, 
                            refdata   = Idents(integrated.mye.seurat),
                            dims      = 1:50)

# Add to the wirka seurat object
mye.wirka.seurat <- AddMetaData(mye.wirka.seurat, metadata = predictions)

# Map to ref UMAP
integrated.mye.seurat <- RunUMAP(integrated.mye.seurat, dims = 1:50, return.model = T)
mye.wirka.seurat <- MapQuery(anchorset       = mye.wirka.anchors, 
                           reference           = integrated.mye.seurat, 
                           query               = mye.wirka.seurat,
                           refdata             = Idents(integrated.mye.seurat), 
                           reference.reduction = "pca", 
                           reduction.model     = "umap")

# Visualize
customUMAP(object = mye.wirka.seurat, group.by = "predicted.id", pt.size = 1.5, label = F, title = "wirka et al. ref mapped", legend.pos = "right",
           shuffle = T, seed = 666, file.name = "wirka_validation_results/wirka.mye_ref_map.UMAP.pdf", cols = M.int_refined.pop.colors, plot.width = 15)
customUMAP(object = mye.wirka.seurat, group.by = "predicted.id", reduction = "ref.umap", pt.size = 1.5, label = F, title = "wirka et al. ref UMAP mapped",
           shuffle = T, seed = 666, file.name = "wirka_validation_results/wirka.mye_ref_umap.UMAP.pdf", legend.pos = "right", cols = M.int_refined.pop.colors, plot.width = 15)
customUMAP(object = integrated.mye.seurat,  pt.size = 1.5, label = F, legend.pos = "right", 
           shuffle = T, seed = 666, file.name = "wirka_validation_results/Mye.UMAP.pdf", cols = M.int_refined.pop.colors, plot.width = 15)
bunchOfCustomPlots(features   = c("TREM2", "LYVE1", "OLR1",  "IL1B",  "CD14", 
                                  "CD9",   "MRC1",  "ABCA1", "TNF",   "FCGR1A",
                                  "GPNMB", "FOLR2", "ABCG1", "CASP1", "FCGR3A"), 
                   object     = mye.wirka.seurat, 
                   group.by   = "predicted.id",
                   Vln.color  = M.int_refined.pop.colors,
                   ncol       = 5,
                   name       = "wirka_validation_results/LAM - resident - foam - inflammatory - general markers pop colors",
                   Vln.draw.names = F,
                   Vln.width  = 30, Vln.height = 15
)

stratifyByExpression(object = mye.wirka.seurat, strat.by = "OLR1",  file.name = "wirka_validation_results/Mye.UMAP.OLR1.pdf",  onlyUMAP = T)
stratifyByExpression(object = mye.wirka.seurat, strat.by = "ABCA1", file.name = "wirka_validation_results/Mye.UMAP.ABCA1.pdf", onlyUMAP = T)
stratifyByExpression(object = mye.wirka.seurat, strat.by = "ABCG1", file.name = "wirka_validation_results/Mye.UMAP.ABCG1.pdf", onlyUMAP = T)
FeaturePlot(object = mye.wirka.seurat, features = c("OLR1", "ABCA1", "ABCG1"), pt.size = 1, order = T)
ggsave("wirka_validation_results/Mye.UMAP.foamy.pdf")
FeaturePlot(object = integrated.mye.seurat, features = c("OLR1", "ABCA1", "ABCG1"), pt.size = 1, order = T)
ggsave("wirka_validation_results/Mye.our_Data.UMAP.foamy.pdf")
