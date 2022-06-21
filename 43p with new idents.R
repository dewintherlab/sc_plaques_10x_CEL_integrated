##==================================================================================
## Use new 10X - CEL-seq integrated populations with the 43p object (mainly for metadata stuff as we don't have it for the 10X guys)
##==================================================================================
dir.create("43p_with_new_refined_idents", showWarnings = F)

# Get new refined idents
int.idents <- Idents(integrated.mye.seurat)
mye.43p.seurat <- AddMetaData(mye.43p.seurat, metadata = int.idents[names(int.idents) %in% colnames(mye.43p.seurat)], col.name = "method.int.idents")

# Plot
customUMAP(object     = mye.43p.seurat, 
           group.by   = "method.int.idents", 
           cols       = M.int_refined.pop.colors, 
           title      = "43p refined populations", 
           legend.pos = "none", 
           plot.width = 15,
           file.name  = "43p_with_new_refined_idents/UMAP 43p integrated populations.pdf" )

customUMAP(object     = mye.43p.seurat, 
           group.by   = "method.int.idents", 
           cols       = M.int_refined.pop.colors, 
           title      = "43p refined populations", 
           legend.pos = "right", 
           plot.width = 15,
           file.name  = "43p_with_new_refined_idents/UMAP 43p integrated populations legend.pdf" )
