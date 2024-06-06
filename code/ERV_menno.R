#=========================================================================================================================
## Plot ERV genes for Menno
##========================================================================================================================

## Setup dir
dir.create("ERV_Menno/top_ERV", showWarnings = F, recursive = T)
dir.create("ERV_Menno/all_TE", showWarnings = F, recursive = T)

## Load gene lists
erv_genes <- scan(file = "ERV_Menno/ERV_genes.txt", what = "list")
TE_genes  <- scan(file = "ERV_Menno/TEname.txt",    what = "list")

## Load TE inclusive Seurat object
TE.seurat <- readRDS(file = "../10X analyses/plaques_and_pbmcs_project/Seurat_Objects/TE.seurat.integrated.RDS")
TE.seurat

## Project our normal UMAP
# Update model
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSeuratObject(final.pop.call.integrated.mye.seurat)
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSCTAssays(final.pop.call.integrated.mye.seurat.updated.sct)
final.pop.call.integrated.mye.seurat.updated.sct <- SCTransform(final.pop.call.integrated.mye.seurat.updated.sct, verbose = T, conserve.memory = T)

# Find anchors
anchors <- FindTransferAnchors(reference            = final.pop.call.integrated.mye.seurat.updated.sct,
                               query                = TE.seurat,
                               normalization.method = "SCT",
                               reference.reduction  = "pca", 
                               recompute.residuals  = T,
                               dims = 1:30
)

# Stash idents
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = Idents(final.pop.call.integrated.mye.seurat.updated.sct), col.name = "ref.pop.idents")
final.pop.call.from_full.integrated.mac.seurat$archetype
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = )

# Map to reference
TE.seurat <- MapQuery(anchorset      = anchors,
                                 query                = TE.seurat,
                                 reference            = final.pop.call.integrated.mye.seurat.updated.sct,
                                 refdata              = list(ref.pop.idents = "ref.pop.idents"),
                                 reference.reduction = "pca", 
                                 reduction.model     = "umap"
)


# ## Extract the macrophages
# # Get 10X plaque mac IDs from our 'normal' seurat object
# macs <- row.names(final.pop.call.integrated.mye.seurat@meta.data)
# macs <- macs[grep("plaque", macs)]
# macs <- gsub("plaque_", "", macs)
# macs <- gsub("-1_[0-9]$", "", macs)
# macs
# 
# # Get the TE IDs
# TE.macs <- row.names(TE.seurat@meta.data)
# 
# # Remove double cells
# x <- list()
# for (i in macs.tmp){
#   x[[i]] <- grep(i, TE.macs)
# }
# 
# x <- x[lengths(x) > 1]
# x
# names(x)
# 
# y <- list()
# for (i in names(x)){
#   y[[i]] <- TE.macs[grep(i, TE.macs)]
# }
# y
# y <- unlist(lapply(y, function(x)x[1]))
# TE.macs <- TE.macs[! TE.macs %in% y]
# TE.macs <- gsub("-[0-9]_[0-9]_[0-9]$", "", TE.macs)
# TE.macs
# 
# # Match the cells
# sum(TE.macs %in% macs)
# macs.to.keep <- row.names(TE.seurat@meta.data)[TE.macs %in% macs]
# 
# # Subset the seurat object
# TE.macs.seurat <- subset(TE.seurat, cells = macs.to.keep)
# TE.macs.seurat
# 
# ## Carry over some metadata
# macs.umap <- final.pop.call.integrated.mye.seurat[['umap']]@cell.embeddings
# macs.umap <- macs.umap[grep("plaque", row.names(macs.umap)),]
# macs2 <- row.names(macs.umap)
# macs2 <- gsub("plaque_", "", macs2)
# macs2 <- gsub("-1_[0-9]$", "", macs2)
# macs2
# 
# TE.cells <- vector()
# for (i in macs2){
#   TE.cells <- c(TE.cells, grep(i, row.names(TE.macs.seurat@meta.data)))
# }
# row.names(macs.umap) <- row.names(TE.macs.seurat@meta.data)[TE.cells]
# 
# TE.macs.seurat[["umap"]] <- CreateDimReducObject(embeddings = macs.umap, assay = "RNA", key = "UMAP_")

## First check the top ERVs
## Make some plots
# All pops
DefaultAssay(final.pop.call.integrated.mye.seurat) <- "RNA"
customUMAP(object = final.pop.call.integrated.mye.seurat, pt.size = 2, cols = M.int_refined.pop.colors, file.name = "ERV_Menno/top_ERV/UMAP_for_reference.pdf", legend.pos = "right", plot.width = 17)
bunchOfCustomPlots(object = final.pop.call.integrated.mye.seurat, features = erv_genes, Vln.draw.names = F, name = "ERV_Menno/top_ERV/ERV_genes", Vln.width = 10, Vln.color = M.int_refined.pop.colors, Vln.height = 10, ncol = 2, feature.pt.size = 1)
customVln(object = final.pop.call.integrated.mye.seurat, features = erv_genes, draw.names = F, name = "ERV_Menno/top_ERV/ERV_genes - violin plot with dots.pdf", width = 10, pt.size = 1, cols = M.int_refined.pop.colors, height = 10, ncol = 2)

# Archetypes
DefaultAssay(final.pop.call.from_full.integrated.mac.seurat) <- "RNA"
customVln(object = final.pop.call.from_full.integrated.mac.seurat, group.by = "archetype", features = erv_genes, draw.names = F, name = "PRC2_Rosalie/PRC2_genes - archetypes - violin plot.pdf", width = 20, pt.size = 0, cols = archetype.colors, height = 30, ncol = 3)
customVln(object = final.pop.call.from_full.integrated.mac.seurat, group.by = "archetype", features = erv_genes, draw.names = F, name = "PRC2_Rosalie/PRC2_genes - archetypes - violin plot with dots.pdf", width = 20, pt.size = 1, cols = archetype.colors, height = 30, ncol = 3)
customDot(object = final.pop.call.from_full.integrated.mac.seurat, group.by = "archetype", features = erv_genes, name = "PRC2_Rosalie/PRC2_genes - archetypes - dot plot.pdf", width = 20, dot.scale = 20, height = 5)
