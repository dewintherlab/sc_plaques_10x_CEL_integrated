#======================================================================
## Analyze solo_TE annotated transposable elements
#======================================================================
## Read data
## TE locus level
##========================================================================================================================
## Plaque samples

# Read matrixes and create seurat objects
plaque.locus_TE.seurat.list <- list()
for(thePatient in c("g003","g005", "g007")){
  cat("Loading data for:", thePatient, "\n")
  theDir <- list.files(path = "raw_data/soloTE/", pattern = thePatient, full.names = T, recursive = F, ignore.case = T, include.dirs = T)
  theDir <- theDir[grep("_locustes", theDir)]
  plaque.locus_TE.seurat.list[[thePatient]] <- Read10X(data.dir = theDir, strip.suffix = T)
  plaque.locus_TE.seurat.list[[thePatient]] <- CreateSeuratObject(counts = plaque.locus_TE.seurat.list[[thePatient]])
}

#=========================================================================================================================
## Normalise the data
dir.create("soloTE.norm_QC/locustes", showWarnings = F, recursive = T)

##==========================================================================================================
## Clean up, normalize, and cluster until we are satisfied
##  BASE: run an unfiltered analysis as base
## Normalize

# Remove empty cells
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)subset(x, subset = nFeature_RNA > 200))

## Add the gene symbols
for(thePatient in c("g003","g005", "g007")){
  plaque.locus_TE.refseq <- row.names(plaque.locus_TE.seurat.list[[thePatient]])
  plaque.locus_TE.refseq <- gsub("^NR-", "NR_", plaque.locus_TE.refseq)
  plaque.locus_TE.refseq <- gsub("^NM-", "NM_", plaque.locus_TE.refseq)

  # Convert to symbols
  plaque.locus_TE.symbol <- AnnotationDbi::select(org.Hs.eg.db, keys = plaque.locus_TE.refseq, keytype = "REFSEQ", columns = "SYMBOL")
  plaque.locus_TE.symbol[is.na(plaque.locus_TE.symbol$SYMBOL),"SYMBOL"] <- plaque.locus_TE.symbol[is.na(plaque.locus_TE.symbol$SYMBOL),"REFSEQ"]
  
  # Reset Symbols to be to rowname compatible
  plaque.locus_TE.symbol$SYMBOL <- make.names(plaque.locus_TE.symbol$SYMBOL, unique = T)
  
  # Add to Seurat object
  plaque.locus_TE.seurat.list[[thePatient]]@assays$RNA@counts@Dimnames[[1]]     <- plaque.locus_TE.symbol$SYMBOL
  plaque.locus_TE.seurat.list[[thePatient]]@assays$RNA@data@Dimnames[[1]]       <- plaque.locus_TE.symbol$SYMBOL
  row.names(plaque.locus_TE.seurat.list[[thePatient]]@assays$RNA@meta.features) <- plaque.locus_TE.symbol$SYMBOL
}

# Normalise
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)UpdateSeuratObject(x))
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)DietSeurat(x))
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)NormalizeData(x, verbose = T))
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)FindVariableFeatures(x, verbose = T))
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)ScaleData(x, verbose = T))

# run sctransform
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)SCTransform(x, verbose = T))

# Save the object
saveRDS(plaque.locus_TE.seurat.list, "Seurat_Objects/main.plaque.locus_TE.seurat.list.prefilter.RDS")
#plaque.locus_TE.seurat.list <- readRDS(file = "Seurat_Objects/main.plaque.locus_TE.seurat.list.prefilter.RDS")

# Run PCA and UMAP embedding
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)RunPCA( x, verbose = T))
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Show some stats
plaque.locus_TE.seurat.list

## Clustering
set.seed(1)
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
plaque.locus_TE.seurat.list <- lapply(plaque.locus_TE.seurat.list, function(x)FindClusters( x, verbose = T, resolution = 1, random.seed = 666))

for(theSample in names(plaque.locus_TE.seurat.list)){
  DimPlot(plaque.locus_TE.seurat.list[[theSample]], label = TRUE, pt.size = 2, label.size = 10) + NoLegend() +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  
  ggsave(paste("soloTE.norm_QC/locustes/", theSample, " Unfiltered Round UMAP clusters.pdf", sep = ""), width = 15, height = 15)

}

# Save the object
saveRDS(plaque.locus_TE.seurat.list, "Seurat_Objects/main.plaque.locus_TE.seurat.list.prefilter.RDS")


#=========================================================================================================================
## Integrate SCTransformed, filtered samples
##========================================================================================================================
dir.create("soloTE.integration_QC/locustes", showWarnings = F, recursive = T)

# Prep for integration
features                    <- SelectIntegrationFeatures(object.list          = plaque.locus_TE.seurat.list, nfeatures = 3000)
plaque.locus_TE.seurat.list <- PrepSCTIntegration(       object.list          = plaque.locus_TE.seurat.list, anchor.features = features)
anchors                     <- FindIntegrationAnchors(   object.list          = plaque.locus_TE.seurat.list, 
                                                         normalization.method = "SCT",
                                                         anchor.features      = features, 
                                                         dims                 = 1:30, 
                                                         reduction            = "rpca", 
                                                         k.anchor             = 20)

# Integrate the samples
plaque.locus_TE.seurat.list <- IntegrateData(anchorset           = anchors, 
                                normalization.method = "SCT", 
                                dims                 = 1:30)
# Show stats
plaque.locus_TE.seurat.list
plaque.locus_TE.seurat <- plaque.locus_TE.seurat.list
rm(plaque.locus_TE.seurat.list)


#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
# Run dimensionality reductions
plaque.locus_TE.seurat <- RunPCA( object = plaque.locus_TE.seurat, verbose = FALSE)
plaque.locus_TE.seurat <- RunUMAP(object = plaque.locus_TE.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
plaque.locus_TE.seurat <- FindNeighbors(object = plaque.locus_TE.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
plaque.locus_TE.seurat <- FindClusters( object = plaque.locus_TE.seurat, resolution = 1.5, verbose = T, random.seed = 666)

# Visualize
DimPlot(plaque.locus_TE.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("soloTE.integration_QC/locustes/Integration cluster UMAP.pdf", width = 10, height = 10)


#=========================================================================================================================
## Isolate myeloid clusters
##========================================================================================================================
dir.create("soloTE.myeloid/locustes", showWarnings = F, recursive = T)

# First check some meyloid genes to see if the (visual) clusters are a fit
mye.gene.symbols <- c("CD68", "CD4", "TREM1", "TREM2", "PLTP", "SELL", "FCGR3A")
bunchOfCustomPlots(object = plaque.locus_TE.seurat, features = mye.gene.symbols, name = "soloTE.integration_QC/locustes/Mye genes")

## Extract macrophages
# Subset myeloid clusters
mye.plaque.locus_TE.seurat <- subset(plaque.locus_TE.seurat, idents = c(6, 9, 10, 12))

# Remove DCs
DefaultAssay(mye.plaque.locus_TE.seurat) <- "RNA"
mye.plaque.locus_TE.seurat <- subset(mye.plaque.locus_TE.seurat, subset = CD1C == 0 & CLEC10A == 0)

#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
mye.plaque.locus_TE.seurat <- SCTransform(mye.plaque.locus_TE.seurat, verbose = T)

# Run dimensionality reductions
mye.plaque.locus_TE.seurat <- RunPCA( object = mye.plaque.locus_TE.seurat, verbose = FALSE)
mye.plaque.locus_TE.seurat <- RunUMAP(object = mye.plaque.locus_TE.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
mye.plaque.locus_TE.seurat <- FindNeighbors(object = mye.plaque.locus_TE.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
mye.plaque.locus_TE.seurat <- FindClusters( object = mye.plaque.locus_TE.seurat, resolution = 1, verbose = T, random.seed = 666)

# Visualize
customUMAP(object = mye.plaque.locus_TE.seurat, file.name = "soloTE.myeloid/locustes/mye UMAP.pdf")

# First check some meyloid genes to see if the (visual) clusters are a fit
bunchOfCustomPlots(object = mye.plaque.locus_TE.seurat, features = mye.gene.symbols, name = "soloTE.myeloid/locustes/Mye genes")


#=========================================================================================================================
## Ref map our idents
##========================================================================================================================
## Project our normal UMAP
DefaultAssay(mye.plaque.locus_TE.seurat) <- "SCT"

# Update model
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSeuratObject(final.pop.call.integrated.mye.seurat)
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSCTAssays(final.pop.call.integrated.mye.seurat.updated.sct)
final.pop.call.integrated.mye.seurat.updated.sct <- SCTransform(final.pop.call.integrated.mye.seurat.updated.sct, verbose = T, conserve.memory = T)
DefaultAssay(final.pop.call.integrated.mye.seurat.updated.sct) <- "SCT"

# Find anchors
anchors <- FindTransferAnchors(reference            = final.pop.call.integrated.mye.seurat.updated.sct,
                               query                = mye.plaque.locus_TE.seurat,
                               normalization.method = "SCT",
                               reference.reduction  = "pca", 
                               recompute.residuals  = T,
                               dims = 1:30, verbose = T
)

# Stash idents
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = Idents(final.pop.call.integrated.mye.seurat.updated.sct), col.name = "ref.pop.idents")
final.pop.call.from_full.integrated.mac.seurat$archetype
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = final.pop.call.from_full.integrated.mac.seurat$archetype, col.name = "archetype")

# Map to reference
mye.plaque.locus_TE.seurat <- MapQuery(anchorset      = anchors,
                                query                 = mye.plaque.locus_TE.seurat,
                                reference             = final.pop.call.integrated.mye.seurat.updated.sct,
                                refdata               = list(ref.pop.idents = "ref.pop.idents", archetype = "archetype"),
                                reference.reduction   = "pca", 
                                reduction.model       = "umap"
)

customUMAP(object = mye.plaque.locus_TE.seurat, reduction = "ref.umap", group.by = "predicted.ref.pop.idents", file.name = "soloTE.myeloid/locustes/mye ref map UMAP.pdf",           legend.pos = "right", cols = M.int_refined.pop.colors, plot.width = 15)
customUMAP(object = mye.plaque.locus_TE.seurat, reduction = "ref.umap", group.by = "predicted.archetype",      file.name = "soloTE.myeloid/locustes/mye ref map archetype UMAP.pdf", legend.pos = "right", cols = archetype.colors)


#=========================================================================================================================
## Call differential TEs
##========================================================================================================================
# Normalise the RNA assay
DefaultAssay(mye.plaque.locus_TE.seurat) <- "RNA"
mye.plaque.locus_TE.seurat <- NormalizeData(mye.plaque.locus_TE.seurat)
mye.plaque.locus_TE.seurat <- FindVariableFeatures(mye.plaque.locus_TE.seurat)
mye.plaque.locus_TE.seurat <- ScaleData(mye.plaque.locus_TE.seurat)

# Find marker TEs
TE.list  <- row.names(mye.plaque.locus_TE.seurat)[grep("Solo", row.names(mye.plaque.locus_TE.seurat))]
ERV.list <- TE.list[grep("ERV", TE.list)]
mye.plaque.locus_TE.markers  <- FindAllMarkers(mye.plaque.locus_TE.seurat, features = TE.list,  only.pos = T)
mye.plaque.locus_ERV.markers <- FindAllMarkers(mye.plaque.locus_TE.seurat, features = ERV.list, only.pos = T)

# Save the top9 markers per cluster
top9.TE.markers  <- mye.plaque.locus_TE.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)
top9.ERV.markers <- mye.plaque.locus_ERV.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# Plot!
for (i in unique(top9.ERV.markers$cluster)){
  bunchOfCustomPlots(mye.plaque.locus_TE.seurat, group.by ="predicted.archetype", name = paste("soloTE.myeloid/locustes/ERV locus cluster ", i, " - top ERV", sep = ""), features = top9.ERV.markers[which(top9.ERV.markers$cluster == i),]$gene, feature.pt.size = 2, Vln.width = 40, Dot.height = 20, Dot.width = 10, Vln.pt.size = 1, Vln.color = archetype.colors, ncol = 3)
}

for(i in unique(top9.TE.markers$cluster)){
  bunchOfCustomPlots(mye.plaque.locus_TE.seurat, group.by ="predicted.archetype", name = paste("soloTE.myeloid/locustes/TE locus cluster ", i, " - top TE", sep = ""), features = top9.TE.markers[which(top9.TE.markers$cluster == i),]$gene, feature.pt.size = 2, Vln.width = 40, Dot.height = 20, Dot.width = 10,  Vln.pt.size = 1, Vln.color = archetype.colors, ncol = 3)
}

#======================================================================
## Read data
## TE class level
##========================================================================================================================
## Plaque samples

# Read matrixes and create seurat objects
plaque.class_TE.seurat.list <- list()
for(thePatient in c("g003","g005", "g007")){
  cat("Loading data for:", thePatient, "\n")
  theDir <- list.files(path = "raw_data/soloTE/", pattern = thePatient, full.names = T, recursive = F, ignore.case = T, include.dirs = T)
  theDir <- theDir[grep("_classtes", theDir)]
  plaque.class_TE.seurat.list[[thePatient]] <- Read10X(data.dir = theDir, strip.suffix = T)
  plaque.class_TE.seurat.list[[thePatient]] <- CreateSeuratObject(counts = plaque.class_TE.seurat.list[[thePatient]])
}

#=========================================================================================================================
## Normalise the data
dir.create("soloTE.norm_QC/classtes", showWarnings = F, recursive = T)

##==========================================================================================================
## Clean up, normalize, and cluster until we are satisfied
##  BASE: run an unfiltered analysis as base
## Normalize

# Remove empty cells
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)subset(x, subset = nFeature_RNA > 200))

## Add the gene symbols
for(thePatient in c("g003","g005", "g007")){
  plaque.class_TE.refseq <- row.names(plaque.class_TE.seurat.list[[thePatient]])
  plaque.class_TE.refseq <- gsub("^NR-", "NR_", plaque.class_TE.refseq)
  plaque.class_TE.refseq <- gsub("^NM-", "NM_", plaque.class_TE.refseq)
  
  # Convert to symbols
  plaque.class_TE.symbol <- AnnotationDbi::select(org.Hs.eg.db, keys = plaque.class_TE.refseq, keytype = "REFSEQ", columns = "SYMBOL")
  plaque.class_TE.symbol[is.na(plaque.class_TE.symbol$SYMBOL),"SYMBOL"] <- plaque.class_TE.symbol[is.na(plaque.class_TE.symbol$SYMBOL),"REFSEQ"]
  
  # Reset Symbols to be to rowname compatible
  plaque.class_TE.symbol$SYMBOL <- make.names(plaque.class_TE.symbol$SYMBOL, unique = T)
  
  # Add to Seurat object
  plaque.class_TE.seurat.list[[thePatient]]@assays$RNA@counts@Dimnames[[1]]     <- plaque.class_TE.symbol$SYMBOL
  plaque.class_TE.seurat.list[[thePatient]]@assays$RNA@data@Dimnames[[1]]       <- plaque.class_TE.symbol$SYMBOL
  row.names(plaque.class_TE.seurat.list[[thePatient]]@assays$RNA@meta.features) <- plaque.class_TE.symbol$SYMBOL
}

# Normalise
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)UpdateSeuratObject(x))
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)DietSeurat(x))
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)NormalizeData(x, verbose = T))
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)FindVariableFeatures(x, verbose = T))
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)ScaleData(x, verbose = T))

# run sctransform
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)SCTransform(x, verbose = T))

# Save the object
saveRDS(plaque.class_TE.seurat.list, "Seurat_Objects/main.plaque.class_TE.seurat.list.prefilter.RDS")
#plaque.class_TE.seurat.list <- readRDS(file = "Seurat_Objects/main.plaque.class_TE.seurat.list.prefilter.RDS")

# Run PCA and UMAP embedding
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)RunPCA( x, verbose = T))
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Show some stats
plaque.class_TE.seurat.list

## Clustering
set.seed(1)
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
plaque.class_TE.seurat.list <- lapply(plaque.class_TE.seurat.list, function(x)FindClusters( x, verbose = T, resolution = 1, random.seed = 666))

for(theSample in names(plaque.class_TE.seurat.list)){
  DimPlot(plaque.class_TE.seurat.list[[theSample]], label = TRUE, pt.size = 2, label.size = 10) + NoLegend() +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  
  ggsave(paste("soloTE.norm_QC/classtes/", theSample, " Unfiltered Round UMAP clusters.pdf", sep = ""), width = 15, height = 15)
  
}

# Save the object
saveRDS(plaque.class_TE.seurat.list, "Seurat_Objects/main.plaque.class_TE.seurat.list.prefilter.RDS")


#=========================================================================================================================
## Integrate SCTransformed, filtered samples
##========================================================================================================================
dir.create("soloTE.integration_QC/classtes", showWarnings = F, recursive = T)

# Prep for integration
features                    <- SelectIntegrationFeatures(object.list          = plaque.class_TE.seurat.list, nfeatures = 3000)
plaque.class_TE.seurat.list <- PrepSCTIntegration(       object.list          = plaque.class_TE.seurat.list, anchor.features = features)
anchors                     <- FindIntegrationAnchors(   object.list          = plaque.class_TE.seurat.list, 
                                                         normalization.method = "SCT",
                                                         anchor.features      = features, 
                                                         dims                 = 1:30, 
                                                         reduction            = "rpca", 
                                                         k.anchor             = 20)

# Integrate the samples
plaque.class_TE.seurat.list <- IntegrateData(anchorset           = anchors, 
                                             normalization.method = "SCT", 
                                             dims                 = 1:30)
# Show stats
plaque.class_TE.seurat.list
plaque.class_TE.seurat <- plaque.class_TE.seurat.list
rm(plaque.class_TE.seurat.list)


#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
# Run dimensionality reductions
plaque.class_TE.seurat <- RunPCA( object = plaque.class_TE.seurat, verbose = FALSE)
plaque.class_TE.seurat <- RunUMAP(object = plaque.class_TE.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
plaque.class_TE.seurat <- FindNeighbors(object = plaque.class_TE.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
plaque.class_TE.seurat <- FindClusters( object = plaque.class_TE.seurat, resolution = 1.5, verbose = T, random.seed = 666)

# Visualize
DimPlot(plaque.class_TE.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("soloTE.integration_QC/classtes/Integration cluster UMAP.pdf", width = 10, height = 10)


#=========================================================================================================================
## Isolate myeloid clusters
##========================================================================================================================
dir.create("soloTE.myeloid/classtes", showWarnings = F, recursive = T)

# First check some meyloid genes to see if the (visual) clusters are a fit
mye.gene.symbols <- c("CD68", "CD4", "TREM1", "TREM2", "PLTP", "SELL", "FCGR3A")
bunchOfCustomPlots(object = plaque.class_TE.seurat, features = mye.gene.symbols, name = "soloTE.integration_QC/classtes/Mye genes")

## Extract macrophages
# Subset myeloid clusters
mye.plaque.class_TE.seurat <- subset(plaque.class_TE.seurat, idents = c(5, 7, 10, 13))

# Remove DCs
DefaultAssay(mye.plaque.class_TE.seurat) <- "RNA"
mye.plaque.class_TE.seurat <- subset(mye.plaque.class_TE.seurat, subset = CD1C == 0 & CLEC10A == 0)

#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
mye.plaque.class_TE.seurat <- SCTransform(mye.plaque.class_TE.seurat, verbose = T)

# Run dimensionality reductions
mye.plaque.class_TE.seurat <- RunPCA( object = mye.plaque.class_TE.seurat, verbose = FALSE)
mye.plaque.class_TE.seurat <- RunUMAP(object = mye.plaque.class_TE.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
mye.plaque.class_TE.seurat <- FindNeighbors(object = mye.plaque.class_TE.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
mye.plaque.class_TE.seurat <- FindClusters( object = mye.plaque.class_TE.seurat, resolution = 1, verbose = T, random.seed = 666)

# Visualize
customUMAP(object = mye.plaque.class_TE.seurat, file.name = "soloTE.myeloid/classtes/mye UMAP.pdf")

# First check some meyloid genes to see if the (visual) clusters are a fit
bunchOfCustomPlots(object = mye.plaque.class_TE.seurat, features = mye.gene.symbols, name = "soloTE.myeloid/classtes/Mye genes")


#=========================================================================================================================
## Ref map our idents
##========================================================================================================================
## Project our normal UMAP
DefaultAssay(mye.plaque.class_TE.seurat) <- "SCT"

# Update model
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSeuratObject(final.pop.call.integrated.mye.seurat)
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSCTAssays(final.pop.call.integrated.mye.seurat.updated.sct)
final.pop.call.integrated.mye.seurat.updated.sct <- SCTransform(final.pop.call.integrated.mye.seurat.updated.sct, verbose = T, conserve.memory = T)
DefaultAssay(final.pop.call.integrated.mye.seurat.updated.sct) <- "SCT"

# Find anchors
anchors <- FindTransferAnchors(reference            = final.pop.call.integrated.mye.seurat.updated.sct,
                               query                = mye.plaque.class_TE.seurat,
                               normalization.method = "SCT",
                               reference.reduction  = "pca", 
                               recompute.residuals  = T,
                               dims = 1:30, verbose = T
)

# Stash idents
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = Idents(final.pop.call.integrated.mye.seurat.updated.sct), col.name = "ref.pop.idents")
final.pop.call.from_full.integrated.mac.seurat$archetype
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = final.pop.call.from_full.integrated.mac.seurat$archetype, col.name = "archetype")

# Map to reference
mye.plaque.class_TE.seurat <- MapQuery(anchorset      = anchors,
                                       query                 = mye.plaque.class_TE.seurat,
                                       reference             = final.pop.call.integrated.mye.seurat.updated.sct,
                                       refdata               = list(ref.pop.idents = "ref.pop.idents", archetype = "archetype"),
                                       reference.reduction   = "pca", 
                                       reduction.model       = "umap"
)

customUMAP(object = mye.plaque.class_TE.seurat, reduction = "ref.umap", group.by = "predicted.ref.pop.idents", file.name = "soloTE.myeloid/classtes/mye ref map UMAP.pdf",           legend.pos = "right", cols = M.int_refined.pop.colors, plot.width = 15)
customUMAP(object = mye.plaque.class_TE.seurat, reduction = "ref.umap", group.by = "predicted.archetype",      file.name = "soloTE.myeloid/classtes/mye ref map archetype UMAP.pdf", legend.pos = "right", cols = archetype.colors)


#=========================================================================================================================
## Call differential TEs
##========================================================================================================================
# Normalise the RNA assay
DefaultAssay(mye.plaque.class_TE.seurat) <- "RNA"
mye.plaque.class_TE.seurat <- NormalizeData(mye.plaque.class_TE.seurat)
mye.plaque.class_TE.seurat <- FindVariableFeatures(mye.plaque.class_TE.seurat)
mye.plaque.class_TE.seurat <- ScaleData(mye.plaque.class_TE.seurat)

# Find marker TEs
TE.list  <- row.names(mye.plaque.class_TE.seurat)[grep("Solo", row.names(mye.plaque.class_TE.seurat))]
#ERV.list <- TE.list[grep("ERV", TE.list)]
#mye.plaque.class_TE.markers  <- FindAllMarkers(mye.plaque.class_TE.seurat, features = TE.list,  only.pos = T)
#mye.plaque.class_ERV.markers <- FindAllMarkers(mye.plaque.class_TE.seurat, features = ERV.list, only.pos = T)

# Save the top9 markers per cluster
#top9.TE.markers <- mye.plaque.class_TE.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# Plot!
bunchOfCustomPlots(mye.plaque.class_TE.seurat, group.by ="predicted.archetype", name = "soloTE.myeloid/classtes/TE classes", features = TE.list, feature.pt.size = 2, Vln.pt.size = 1, Vln.color = archetype.colors, ncol = 3)



#======================================================================
## Read data
## TE family level
##========================================================================================================================
## Plaque samples

# Read matrixes and create seurat objects
plaque.family_TE.seurat.list <- list()
for(thePatient in c("g003","g005", "g007")){
  cat("Loading data for:", thePatient, "\n")
  theDir <- list.files(path = "raw_data/soloTE/", pattern = thePatient, full.names = T, recursive = F, ignore.case = T, include.dirs = T)
  theDir <- theDir[grep("_familytes", theDir)]
  plaque.family_TE.seurat.list[[thePatient]] <- Read10X(data.dir = theDir, strip.suffix = T)
  plaque.family_TE.seurat.list[[thePatient]] <- CreateSeuratObject(counts = plaque.family_TE.seurat.list[[thePatient]])
}

#=========================================================================================================================
## Normalise the data
dir.create("soloTE.norm_QC/familytes",   showWarnings = F, recursive = T)

##==========================================================================================================
## Clean up, normalize, and cluster until we are satisfied
##  BASE: run an unfiltered analysis as base
## Normalize

# Remove empty cells
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)subset(x, subset = nFeature_RNA > 200))

## Add the gene symbols
for(thePatient in c("g003","g005", "g007")){
  plaque.family_TE.refseq <- row.names(plaque.family_TE.seurat.list[[thePatient]])
  plaque.family_TE.refseq <- gsub("^NR-", "NR_", plaque.family_TE.refseq)
  plaque.family_TE.refseq <- gsub("^NM-", "NM_", plaque.family_TE.refseq)
  
  # Convert to symbols
  plaque.family_TE.symbol <- AnnotationDbi::select(org.Hs.eg.db, keys = plaque.family_TE.refseq, keytype = "REFSEQ", columns = "SYMBOL")
  plaque.family_TE.symbol[is.na(plaque.family_TE.symbol$SYMBOL),"SYMBOL"] <- plaque.family_TE.symbol[is.na(plaque.family_TE.symbol$SYMBOL),"REFSEQ"]
  
  # Reset Symbols to be to rowname compatible
  plaque.family_TE.symbol$SYMBOL <- make.names(plaque.family_TE.symbol$SYMBOL, unique = T)
  
  # Add to Seurat object
  plaque.family_TE.seurat.list[[thePatient]]@assays$RNA@counts@Dimnames[[1]]     <- plaque.family_TE.symbol$SYMBOL
  plaque.family_TE.seurat.list[[thePatient]]@assays$RNA@data@Dimnames[[1]]       <- plaque.family_TE.symbol$SYMBOL
  row.names(plaque.family_TE.seurat.list[[thePatient]]@assays$RNA@meta.features) <- plaque.family_TE.symbol$SYMBOL
}

# Normalise
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)UpdateSeuratObject(x))
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)DietSeurat(x))
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)NormalizeData(x, verbose = T))
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)FindVariableFeatures(x, verbose = T))
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)ScaleData(x, verbose = T))

# run sctransform
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)SCTransform(x, verbose = T))

# Save the object
saveRDS(plaque.family_TE.seurat.list, "Seurat_Objects/main.plaque.family_TE.seurat.list.prefilter.RDS")
#plaque.family_TE.seurat.list <- readRDS(file = "Seurat_Objects/main.plaque.family_TE.seurat.list.prefilter.RDS")

# Run PCA and UMAP embedding
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)RunPCA( x, verbose = T))
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Show some stats
plaque.family_TE.seurat.list

## Clustering
set.seed(1)
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
plaque.family_TE.seurat.list <- lapply(plaque.family_TE.seurat.list, function(x)FindClusters( x, verbose = T, resolution = 1, random.seed = 666))

for(theSample in names(plaque.family_TE.seurat.list)){
  DimPlot(plaque.family_TE.seurat.list[[theSample]], label = TRUE, pt.size = 2, label.size = 10) + NoLegend() +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  
  ggsave(paste("soloTE.norm_QC/familytes/", theSample, " Unfiltered Round UMAP clusters.pdf", sep = ""), width = 15, height = 15)
  
}

# Save the object
saveRDS(plaque.family_TE.seurat.list, "Seurat_Objects/main.plaque.family_TE.seurat.list.prefilter.RDS")


#=========================================================================================================================
## Integrate SCTransformed, filtered samples
##========================================================================================================================
dir.create("soloTE.integration_QC/familytes", showWarnings = F, recursive = T)

# Prep for integration
features                    <- SelectIntegrationFeatures(object.list          = plaque.family_TE.seurat.list, nfeatures = 3000)
plaque.family_TE.seurat.list <- PrepSCTIntegration(       object.list          = plaque.family_TE.seurat.list, anchor.features = features)
anchors                     <- FindIntegrationAnchors(   object.list          = plaque.family_TE.seurat.list, 
                                                         normalization.method = "SCT",
                                                         anchor.features      = features, 
                                                         dims                 = 1:30, 
                                                         reduction            = "rpca", 
                                                         k.anchor             = 20)

# Integrate the samples
plaque.family_TE.seurat.list <- IntegrateData(anchorset           = anchors, 
                                             normalization.method = "SCT", 
                                             dims                 = 1:30)
# Show stats
plaque.family_TE.seurat.list
plaque.family_TE.seurat <- plaque.family_TE.seurat.list
rm(plaque.family_TE.seurat.list)


#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
# Run dimensionality reductions
plaque.family_TE.seurat <- RunPCA( object = plaque.family_TE.seurat, verbose = FALSE)
plaque.family_TE.seurat <- RunUMAP(object = plaque.family_TE.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
plaque.family_TE.seurat <- FindNeighbors(object = plaque.family_TE.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
plaque.family_TE.seurat <- FindClusters( object = plaque.family_TE.seurat, resolution = 1.5, verbose = T, random.seed = 666)

# Visualize
DimPlot(plaque.family_TE.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("soloTE.integration_QC/familytes/Integration cluster UMAP.pdf", width = 10, height = 10)


#=========================================================================================================================
## Isolate myeloid clusters
##========================================================================================================================
dir.create("soloTE.myeloid/familytes", showWarnings = F, recursive = T)

# First check some meyloid genes to see if the (visual) clusters are a fit
mye.gene.symbols <- c("CD68", "CD4", "TREM1", "TREM2", "PLTP", "SELL", "FCGR3A")
bunchOfCustomPlots(object = plaque.family_TE.seurat, features = mye.gene.symbols, name = "soloTE.integration_QC/familytes/Mye genes")

## Extract macrophages
# Subset myeloid clusters
mye.plaque.family_TE.seurat <- subset(plaque.family_TE.seurat, idents = c(4, 7, 11, 13))

# Remove DCs
DefaultAssay(mye.plaque.family_TE.seurat) <- "RNA"
mye.plaque.family_TE.seurat <- subset(mye.plaque.family_TE.seurat, subset = CD1C == 0 & CLEC10A == 0)

#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
mye.plaque.family_TE.seurat <- SCTransform(mye.plaque.family_TE.seurat, verbose = T)

# Run dimensionality reductions
mye.plaque.family_TE.seurat <- RunPCA( object = mye.plaque.family_TE.seurat, verbose = FALSE)
mye.plaque.family_TE.seurat <- RunUMAP(object = mye.plaque.family_TE.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
mye.plaque.family_TE.seurat <- FindNeighbors(object = mye.plaque.family_TE.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
mye.plaque.family_TE.seurat <- FindClusters( object = mye.plaque.family_TE.seurat, resolution = 1, verbose = T, random.seed = 666)

# Visualize
customUMAP(object = mye.plaque.family_TE.seurat, file.name = "soloTE.myeloid/familytes/mye UMAP.pdf")

# First check some meyloid genes to see if the (visual) clusters are a fit
bunchOfCustomPlots(object = mye.plaque.family_TE.seurat, features = mye.gene.symbols, name = "soloTE.myeloid/familytes/Mye genes")


#=========================================================================================================================
## Ref map our idents
##========================================================================================================================
## Project our normal UMAP
DefaultAssay(mye.plaque.family_TE.seurat) <- "SCT"

# Update model
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSeuratObject(final.pop.call.integrated.mye.seurat)
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSCTAssays(final.pop.call.integrated.mye.seurat.updated.sct)
final.pop.call.integrated.mye.seurat.updated.sct <- SCTransform(final.pop.call.integrated.mye.seurat.updated.sct, verbose = T, conserve.memory = T)
DefaultAssay(final.pop.call.integrated.mye.seurat.updated.sct) <- "SCT"

# Find anchors
anchors <- FindTransferAnchors(reference            = final.pop.call.integrated.mye.seurat.updated.sct,
                               query                = mye.plaque.family_TE.seurat,
                               normalization.method = "SCT",
                               reference.reduction  = "pca", 
                               recompute.residuals  = T,
                               dims = 1:30, verbose = T
)

# Stash idents
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = Idents(final.pop.call.integrated.mye.seurat.updated.sct), col.name = "ref.pop.idents")
final.pop.call.from_full.integrated.mac.seurat$archetype
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = final.pop.call.from_full.integrated.mac.seurat$archetype, col.name = "archetype")

# Map to reference
mye.plaque.family_TE.seurat <- MapQuery(anchorset      = anchors,
                                       query                 = mye.plaque.family_TE.seurat,
                                       reference             = final.pop.call.integrated.mye.seurat.updated.sct,
                                       refdata               = list(ref.pop.idents = "ref.pop.idents", archetype = "archetype"),
                                       reference.reduction   = "pca", 
                                       reduction.model       = "umap"
)

customUMAP(object = mye.plaque.family_TE.seurat, reduction = "ref.umap", group.by = "predicted.ref.pop.idents", file.name = "soloTE.myeloid/familytes/mye ref map UMAP.pdf",           legend.pos = "right", cols = M.int_refined.pop.colors, plot.width = 15)
customUMAP(object = mye.plaque.family_TE.seurat, reduction = "ref.umap", group.by = "predicted.archetype",      file.name = "soloTE.myeloid/familytes/mye ref map archetype UMAP.pdf", legend.pos = "right", cols = archetype.colors)


#=========================================================================================================================
## Call differential TEs
##========================================================================================================================
# Normalise the RNA assay
DefaultAssay(mye.plaque.family_TE.seurat) <- "RNA"
mye.plaque.family_TE.seurat <- NormalizeData(mye.plaque.family_TE.seurat)
mye.plaque.family_TE.seurat <- FindVariableFeatures(mye.plaque.family_TE.seurat)
mye.plaque.family_TE.seurat <- ScaleData(mye.plaque.family_TE.seurat)

# Find marker TEs
TE.list  <- row.names(mye.plaque.family_TE.seurat)[grep("Solo", row.names(mye.plaque.family_TE.seurat))]
ERV.list <- TE.list[grep("ERV", TE.list)]
mye.plaque.family_TE.markers  <- FindAllMarkers(mye.plaque.family_TE.seurat, features = TE.list,  only.pos = T)
mye.plaque.family_ERV.markers <- FindAllMarkers(mye.plaque.family_TE.seurat, features = ERV.list, only.pos = T)

# Save the top9 markers per cluster
top9.TE.markers <- mye.plaque.family_TE.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# Plot!
bunchOfCustomPlots(mye.plaque.family_TE.seurat, group.by ="predicted.archetype", name = "soloTE.myeloid/familytes/ERV families", features = ERV.list, feature.pt.size = 2, Vln.pt.size = 1, Vln.color = archetype.colors, ncol = 3)

for(i in unique(top9.TE.markers$cluster)){
  bunchOfCustomPlots(mye.plaque.family_TE.seurat, group.by ="predicted.archetype", name = paste("soloTE.myeloid/familytes/TE familytes cluster ", i, " - top TE", sep = ""), features = top9.TE.markers[which(top9.TE.markers$cluster == i),]$gene, feature.pt.size = 2, Vln.pt.size = 1, Vln.color = archetype.colors, ncol = 3)
}





#======================================================================
## Read data
## TE subfamily level
##========================================================================================================================
## Plaque samples

# Read matrixes and create seurat objects
plaque.subfamily_TE.seurat.list <- list()
for(thePatient in c("g003","g005", "g007")){
  cat("Loading data for:", thePatient, "\n")
  theDir <- list.files(path = "raw_data/soloTE/", pattern = thePatient, full.names = T, recursive = F, ignore.case = T, include.dirs = T)
  theDir <- theDir[grep("_subfamilytes", theDir)]
  plaque.subfamily_TE.seurat.list[[thePatient]] <- Read10X(data.dir = theDir, strip.suffix = T)
  plaque.subfamily_TE.seurat.list[[thePatient]] <- CreateSeuratObject(counts = plaque.subfamily_TE.seurat.list[[thePatient]])
}

#=========================================================================================================================
## Normalise the data
dir.create("soloTE.norm_QC/subfamilytes",   showWarnings = F, recursive = T)

##==========================================================================================================
## Clean up, normalize, and cluster until we are satisfied
##  BASE: run an unfiltered analysis as base
## Normalize

# Remove empty cells
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)subset(x, subset = nFeature_RNA > 200))

## Add the gene symbols
for(thePatient in c("g003","g005", "g007")){
  plaque.subfamily_TE.refseq <- row.names(plaque.subfamily_TE.seurat.list[[thePatient]])
  plaque.subfamily_TE.refseq <- gsub("^NR-", "NR_", plaque.subfamily_TE.refseq)
  plaque.subfamily_TE.refseq <- gsub("^NM-", "NM_", plaque.subfamily_TE.refseq)
  
  # Convert to symbols
  plaque.subfamily_TE.symbol <- AnnotationDbi::select(org.Hs.eg.db, keys = plaque.subfamily_TE.refseq, keytype = "REFSEQ", columns = "SYMBOL")
  plaque.subfamily_TE.symbol[is.na(plaque.subfamily_TE.symbol$SYMBOL),"SYMBOL"] <- plaque.subfamily_TE.symbol[is.na(plaque.subfamily_TE.symbol$SYMBOL),"REFSEQ"]
  
  # Reset Symbols to be to rowname compatible
  plaque.subfamily_TE.symbol$SYMBOL <- make.names(plaque.subfamily_TE.symbol$SYMBOL, unique = T)
  
  # Add to Seurat object
  plaque.subfamily_TE.seurat.list[[thePatient]]@assays$RNA@counts@Dimnames[[1]]     <- plaque.subfamily_TE.symbol$SYMBOL
  plaque.subfamily_TE.seurat.list[[thePatient]]@assays$RNA@data@Dimnames[[1]]       <- plaque.subfamily_TE.symbol$SYMBOL
  row.names(plaque.subfamily_TE.seurat.list[[thePatient]]@assays$RNA@meta.features) <- plaque.subfamily_TE.symbol$SYMBOL
}

# Normalise
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)UpdateSeuratObject(x))
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)DietSeurat(x))
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)NormalizeData(x, verbose = T))
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)FindVariableFeatures(x, verbose = T))
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)ScaleData(x, verbose = T))

# run sctransform
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)SCTransform(x, verbose = T))

# Save the object
saveRDS(plaque.subfamily_TE.seurat.list, "Seurat_Objects/main.plaque.subfamily_TE.seurat.list.prefilter.RDS")
#plaque.subfamily_TE.seurat.list <- readRDS(file = "Seurat_Objects/main.plaque.subfamily_TE.seurat.list.prefilter.RDS")

# Run PCA and UMAP embedding
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)RunPCA( x, verbose = T))
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Show some stats
plaque.subfamily_TE.seurat.list

## Clustering
set.seed(1)
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
plaque.subfamily_TE.seurat.list <- lapply(plaque.subfamily_TE.seurat.list, function(x)FindClusters( x, verbose = T, resolution = 1, random.seed = 666))

for(theSample in names(plaque.subfamily_TE.seurat.list)){
  DimPlot(plaque.subfamily_TE.seurat.list[[theSample]], label = TRUE, pt.size = 2, label.size = 10) + NoLegend() +
    theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
    theme(panel.background = element_blank(),
          title            = element_text(size = 16, face = "bold")
    )
  
  ggsave(paste("soloTE.norm_QC/subfamilytes/", theSample, " Unfiltered Round UMAP clusters.pdf", sep = ""), width = 15, height = 15)
  
}

# Save the object
saveRDS(plaque.subfamily_TE.seurat.list, "Seurat_Objects/main.plaque.subfamily_TE.seurat.list.prefilter.RDS")


#=========================================================================================================================
## Integrate SCTransformed, filtered samples
##========================================================================================================================
dir.create("soloTE.integration_QC/subfamilytes", showWarnings = F, recursive = T)

# Prep for integration
features                    <- SelectIntegrationFeatures(object.list          = plaque.subfamily_TE.seurat.list, nfeatures = 3000)
plaque.subfamily_TE.seurat.list <- PrepSCTIntegration(       object.list          = plaque.subfamily_TE.seurat.list, anchor.features = features)
anchors                     <- FindIntegrationAnchors(   object.list          = plaque.subfamily_TE.seurat.list, 
                                                         normalization.method = "SCT",
                                                         anchor.features      = features, 
                                                         dims                 = 1:30, 
                                                         reduction            = "rpca", 
                                                         k.anchor             = 20)

# Integrate the samples
plaque.subfamily_TE.seurat.list <- IntegrateData(anchorset           = anchors, 
                                              normalization.method = "SCT", 
                                              dims                 = 1:30)
# Show stats
plaque.subfamily_TE.seurat.list
plaque.subfamily_TE.seurat <- plaque.subfamily_TE.seurat.list
rm(plaque.subfamily_TE.seurat.list)


#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
# Run dimensionality reductions
plaque.subfamily_TE.seurat <- RunPCA( object = plaque.subfamily_TE.seurat, verbose = FALSE)
plaque.subfamily_TE.seurat <- RunUMAP(object = plaque.subfamily_TE.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
plaque.subfamily_TE.seurat <- FindNeighbors(object = plaque.subfamily_TE.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
plaque.subfamily_TE.seurat <- FindClusters( object = plaque.subfamily_TE.seurat, resolution = 1.5, verbose = T, random.seed = 666)

# Visualize
DimPlot(plaque.subfamily_TE.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T) +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("soloTE.integration_QC/subfamilytes/Integration cluster UMAP.pdf", width = 10, height = 10)


#=========================================================================================================================
## Isolate myeloid clusters
##========================================================================================================================
dir.create("soloTE.myeloid/subfamilytes", showWarnings = F, recursive = T)

# First check some meyloid genes to see if the (visual) clusters are a fit
mye.gene.symbols <- c("CD68", "CD4", "TREM1", "TREM2", "PLTP", "SELL", "FCGR3A")
bunchOfCustomPlots(object = plaque.subfamily_TE.seurat, features = mye.gene.symbols, name = "soloTE.integration_QC/subfamilytes/Mye genes")

## Extract macrophages
# Subset myeloid clusters
mye.plaque.subfamily_TE.seurat <- subset(plaque.subfamily_TE.seurat, idents = c(4, 7, 9))

# Remove DCs
DefaultAssay(mye.plaque.subfamily_TE.seurat) <- "RNA"
mye.plaque.subfamily_TE.seurat <- subset(mye.plaque.subfamily_TE.seurat, subset = CD1C == 0 & CLEC10A == 0)

#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
mye.plaque.subfamily_TE.seurat <- SCTransform(mye.plaque.subfamily_TE.seurat, verbose = T)

# Run dimensionality reductions
mye.plaque.subfamily_TE.seurat <- RunPCA( object = mye.plaque.subfamily_TE.seurat, verbose = FALSE)
mye.plaque.subfamily_TE.seurat <- RunUMAP(object = mye.plaque.subfamily_TE.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
mye.plaque.subfamily_TE.seurat <- FindNeighbors(object = mye.plaque.subfamily_TE.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
mye.plaque.subfamily_TE.seurat <- FindClusters( object = mye.plaque.subfamily_TE.seurat, resolution = 1, verbose = T, random.seed = 666)

# Visualize
customUMAP(object = mye.plaque.subfamily_TE.seurat, file.name = "soloTE.myeloid/subfamilytes/mye UMAP.pdf")

# First check some meyloid genes to see if the (visual) clusters are a fit
bunchOfCustomPlots(object = mye.plaque.subfamily_TE.seurat, features = mye.gene.symbols, name = "soloTE.myeloid/subfamilytes/Mye genes")


#=========================================================================================================================
## Ref map our idents
##========================================================================================================================
## Project our normal UMAP
DefaultAssay(mye.plaque.subfamily_TE.seurat) <- "SCT"

# Update model
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSeuratObject(final.pop.call.integrated.mye.seurat)
final.pop.call.integrated.mye.seurat.updated.sct <- UpdateSCTAssays(final.pop.call.integrated.mye.seurat.updated.sct)
final.pop.call.integrated.mye.seurat.updated.sct <- SCTransform(final.pop.call.integrated.mye.seurat.updated.sct, verbose = T, conserve.memory = T)
DefaultAssay(final.pop.call.integrated.mye.seurat.updated.sct) <- "SCT"

# Find anchors
anchors <- FindTransferAnchors(reference            = final.pop.call.integrated.mye.seurat.updated.sct,
                               query                = mye.plaque.subfamily_TE.seurat,
                               normalization.method = "SCT",
                               reference.reduction  = "pca", 
                               recompute.residuals  = T,
                               dims = 1:30, verbose = T
)

# Stash idents
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = Idents(final.pop.call.integrated.mye.seurat.updated.sct), col.name = "ref.pop.idents")
final.pop.call.from_full.integrated.mac.seurat$archetype
final.pop.call.integrated.mye.seurat.updated.sct <- AddMetaData(final.pop.call.integrated.mye.seurat.updated.sct, metadata = final.pop.call.from_full.integrated.mac.seurat$archetype, col.name = "archetype")

# Map to reference
mye.plaque.subfamily_TE.seurat <- MapQuery(anchorset      = anchors,
                                        query                 = mye.plaque.subfamily_TE.seurat,
                                        reference             = final.pop.call.integrated.mye.seurat.updated.sct,
                                        refdata               = list(ref.pop.idents = "ref.pop.idents", archetype = "archetype"),
                                        reference.reduction   = "pca", 
                                        reduction.model       = "umap"
)

customUMAP(object = mye.plaque.subfamily_TE.seurat, reduction = "ref.umap", group.by = "predicted.ref.pop.idents", file.name = "soloTE.myeloid/subfamilytes/mye ref map UMAP.pdf",           legend.pos = "right", cols = M.int_refined.pop.colors, plot.width = 15)
customUMAP(object = mye.plaque.subfamily_TE.seurat, reduction = "ref.umap", group.by = "predicted.archetype",      file.name = "soloTE.myeloid/subfamilytes/mye ref map archetype UMAP.pdf", legend.pos = "right", cols = archetype.colors)


#=========================================================================================================================
## Call differential TEs
##========================================================================================================================
# Normalise the RNA assay
DefaultAssay(mye.plaque.subfamily_TE.seurat) <- "RNA"
mye.plaque.subfamily_TE.seurat <- NormalizeData(mye.plaque.subfamily_TE.seurat)
mye.plaque.subfamily_TE.seurat <- FindVariableFeatures(mye.plaque.subfamily_TE.seurat)
mye.plaque.subfamily_TE.seurat <- ScaleData(mye.plaque.subfamily_TE.seurat)

# Find marker TEs
TE.list  <- row.names(mye.plaque.subfamily_TE.seurat)[grep("Solo", row.names(mye.plaque.subfamily_TE.seurat))]
ERV.list <- TE.list[grep("ERV", TE.list)]
mye.plaque.subfamily_TE.markers  <- FindAllMarkers(mye.plaque.subfamily_TE.seurat, features = TE.list,  only.pos = T)
mye.plaque.subfamily_ERV.markers <- FindAllMarkers(mye.plaque.subfamily_TE.seurat, features = ERV.list, only.pos = T)

# Save the top9 markers per cluster
top9.TE.markers <- mye.plaque.subfamily_TE.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# Plot!
length(ERV.list)
for (i in seq(1,54,6)){
  bunchOfCustomPlots(mye.plaque.subfamily_TE.seurat, group.by ="predicted.archetype", name = paste("soloTE.myeloid/subfamilytes/ERV subfamilytes plot number ", i, " through ", i+5, "- not DE", sep = ""), features = ERV.list[seq(i,i+5)], feature.pt.size = 2, Vln.width = 20, Vln.pt.size = 1, Vln.color = archetype.colors, ncol = 3)
}


for(i in unique(top9.TE.markers$cluster)){
  bunchOfCustomPlots(mye.plaque.subfamily_TE.seurat, group.by ="predicted.archetype", name = paste("soloTE.myeloid/subfamilytes/TE subfamilytes cluster ", i, " - top TE", sep = ""), features = top9.TE.markers[which(top9.TE.markers$cluster == i),]$gene, feature.pt.size = 2, Vln.width = 20, Vln.pt.size = 1, Vln.color = archetype.colors, ncol = 3)
}





