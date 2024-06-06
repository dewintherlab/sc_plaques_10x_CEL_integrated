#=========================================================================================================================
## Plot PRC2 genes for Rosalie
##========================================================================================================================

## Setup dir
dir.create("PRC2_Rosalie", showWarnings = F)

## Load gene list
prc2_genes <- scan(file = "PRC2_Rosalie/PRC2_genes.txt", what = "list")

## Make some plots
# All pops
DefaultAssay(final.pop.call.integrated.mye.seurat) <- "RNA"
customUMAP(object = final.pop.call.integrated.mye.seurat, pt.size = 2, cols = M.int_refined.pop.colors, file.name = "PRC2_Rosalie/UMAP_for_reference.pdf", legend.pos = "right", plot.width = 17)
bunchOfCustomPlots(object = final.pop.call.integrated.mye.seurat, features = prc2_genes, Vln.draw.names = F, name = "PRC2_Rosalie/PRC2_genes", Vln.width = 20, Vln.color = M.int_refined.pop.colors, Vln.height = 30, ncol = 3, feature.pt.size = 1)
customVln(object = final.pop.call.integrated.mye.seurat, features = prc2_genes, draw.names = F, name = "PRC2_Rosalie/PRC2_genes - violin plot with dots.pdf", width = 20, pt.size = 1, cols = M.int_refined.pop.colors, height = 30, ncol = 3)

# Archetypes
DefaultAssay(final.pop.call.from_full.integrated.mac.seurat) <- "RNA"
customVln(object = final.pop.call.from_full.integrated.mac.seurat, group.by = "archetype", features = prc2_genes, draw.names = F, name = "PRC2_Rosalie/PRC2_genes - archetypes - violin plot.pdf", width = 20, pt.size = 0, cols = archetype.colors, height = 30, ncol = 3)
customVln(object = final.pop.call.from_full.integrated.mac.seurat, group.by = "archetype", features = prc2_genes, draw.names = F, name = "PRC2_Rosalie/PRC2_genes - archetypes - violin plot with dots.pdf", width = 20, pt.size = 1, cols = archetype.colors, height = 30, ncol = 3)
customDot(object = final.pop.call.from_full.integrated.mac.seurat, group.by = "archetype", features = prc2_genes, name = "PRC2_Rosalie/PRC2_genes - archetypes - dot plot.pdf", width = 20, dot.scale = 20, height = 5)
