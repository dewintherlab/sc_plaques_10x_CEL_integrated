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

## Subset the macrophages
# Retrieve mac cell ids
macs <- row.names(final.pop.call.integrated.mye.velocyto.seurat@meta.data)
macs <- macs[grep("plaque", macs)]
macs <- gsub("plaque_", "", macs)

# Retrieve TE cell ids and normalize
TE.cells <- row.names(TE.seurat@meta.data)
TE.cells <- gsub("_[0-9]$", "", TE.cells)
sum(TE.cells %in% macs)

# Keep only the mac TE.cells
macs.to.keep   <- row.names(TE.seurat@meta.data)[TE.cells %in% macs]
TE.macs.seurat <- subset(TE.seurat, cells = macs.to.keep)

## Move the mac idents over to the new TE object
# Need to normalize cell ids. Macs have 2 digit sufiix, TE cells 3 digits
head(Idents(final.pop.call.integrated.mye.velocyto.seurat))
head(row.names(TE.macs.seurat@meta.data))

# Get mac barcodes without suffix
s.i <- gsub("plaque_", "", names(Idents(final.pop.call.integrated.mye.velocyto.seurat)))
s.i <- gsub("-[0-9]_[0-9]$", "", s.i)

# Compare the barcodes without suffix to see if there is a pattern we can use
for(i in s.i){
  cat(i, row.names(TE.macs.seurat@meta.data)[grep(i, row.names(TE.macs.seurat@meta.data))], "\n")
}

# Turns out all TE.cells with a mathcing barcode have a _2 as the third suffix digits (and a couple double its with a 1), so let's normalise the names by adding a _2
mac.idents <- Idents(final.pop.call.integrated.mye.velocyto.seurat)
names(mac.idents) <- paste0(gsub("plaque_", "", names(Idents(final.pop.call.integrated.mye.velocyto.seurat))), "_2")

# Add to the TE.object
TE.macs.seurat <- AddMetaData(TE.macs.seurat, mac.idents, col.name = "mac.idents")

# Visualise
TE.macs.seurat <- RunPCA(TE.macs.seurat)
TE.macs.seurat <- RunUMAP(TE.macs.seurat, dims = 1:30)
customUMAP(object = TE.macs.seurat, group.by = "mac.idents", pt.size = 3, cols = M.int_refined.pop.colors, file.name = "ERV_Menno/top_ERV/TE.UMAP_space.pdf", legend.pos = "right", plot.width = 15)


## Carry over the original macs UMAP
macs.umap <- final.pop.call.integrated.mye.velocyto.seurat[['umap']]@cell.embeddings
macs.umap <- macs.umap[grep("plaque", row.names(macs.umap)),]

row.names(macs.umap) <- paste0(gsub("plaque_", "", row.names(macs.umap)), "_2")
TE.macs.seurat[["umap"]] <- CreateDimReducObject(embeddings = macs.umap, assay = "RNA", key = "UMAP_")

# Clean up the 8 superfluous cells from the conversion
TE.macs.seurat <- subset(TE.macs.seurat, cells = row.names(TE.macs.seurat@meta.data)[is.na(TE.macs.seurat@meta.data$mac.idents)], invert = T)

# Visualise
customUMAP(object = TE.macs.seurat, group.by = "mac.idents", pt.size = 3, cols = M.int_refined.pop.colors, file.name = "ERV_Menno/top_ERV/Mac.UMAP_space.pdf", legend.pos = "right", plot.width = 15)

## Add archetypes
tmp.idents <- as.vector(TE.macs.seurat$mac.idents)
unique(tmp.idents)
tmp.idents[grep("Resident-like Macrophages", tmp.idents)] <- "Resident-like"
tmp.idents[grep("Lipid Associated Macrophages", tmp.idents)] <- "LAM"
tmp.idents[grep("Foamy", tmp.idents)] <- "iLAM"
tmp.idents[grep("Inflammatory", tmp.idents)] <- "Inflammatory"
TE.macs.seurat <- AddMetaData(TE.macs.seurat, tmp.idents, col.name = "archetype")
customUMAP(object = TE.macs.seurat, group.by = "archetype", pt.size = 3, cols = archetype.colors, file.name = "ERV_Menno/top_ERV/Mac.archetype.UMAP_space.pdf", legend.pos = "right", plot.width = 10)

## Plot ERVs
DefaultAssay(TE.macs.seurat) <- "RNA"

extra.erv.genes <- c("MER101B",
"MLT1B-int",
"MER92B",
"MER67C",
"MLT1F1-int",
"LTR49-int",
"MER57F")

for(theERV in extra.erv.genes){
  cat("Plotting ", theERV, "...\n", sep = "")
  customVln(object = TE.macs.seurat, features = theERV, group.by = "mac.idents", draw.names = T, name = paste("ERV_Menno/top_ERV/", theERV, " - population violin plot with dots.pdf",    sep = ""), , pt.size = 2, cols = M.int_refined.pop.colors, width = 10, height = 12)
  customVln(object = TE.macs.seurat, features = theERV, group.by = "mac.idents", draw.names = T, name = paste("ERV_Menno/top_ERV/", theERV, " - population violin plot without dots.pdf", sep = ""), , pt.size = 0, cols = M.int_refined.pop.colors, width = 10, height = 12)
  customVln(object = TE.macs.seurat, features = theERV, group.by = "archetype", draw.names = T, name = paste("ERV_Menno/top_ERV/", theERV, " - archetype violin plot with dots.pdf",    sep = ""), , pt.size = 2, cols = archetype.colors, width = 10, height = 8)
  customVln(object = TE.macs.seurat, features = theERV, group.by = "archetype", draw.names = T, name = paste("ERV_Menno/top_ERV/", theERV, " - archetype violin plot without dots.pdf", sep = ""), , pt.size = 0, cols = archetype.colors, width = 10, height = 8)
  customFeature(object = TE.macs.seurat, features = theERV, name = paste("ERV_Menno/top_ERV/", theERV, " - feature plot.pdf", sep = ""), pt.size = 3, width = 5, height = 5)
}

