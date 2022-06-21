#=================================================================================================================================
## Compare monocytes and macrophages to PBMC reference data.
## Use ar recent study (GSE168732) that did 10X runs on 3 healthy controls
##================================================================================================================================
dir.create("compare_monocytes_to_PBMC_ref", showWarnings = F)

# Fetch sample dirs
all.sample.dirs <- list.files(path = "raw_data/Wang_GSE168732", include.dirs = T, full.names = T)

# Load raw data
ctrl.samples.raw <- list()
for(theSample in all.sample.dirs){
  cat("Loading CTRL sample:", basename(theSample), "...\n")
  
  # Load from disk
  ctrl.samples.raw[[basename(theSample)]] <- Read10X(data.dir = theSample)
}

# Turn into seurat objects
ctrl.samples.seurat <- list()
for(theSample in names(ctrl.samples.raw)){
  cat("Processing:", theSample, "\n")
  ctrl.samples.seurat[[theSample]] <- CreateSeuratObject(counts       = ctrl.samples.raw[[theSample]],
                                                         project      = "CTRL", 
                                                         min.cells    = 3,
                                                         min.features = 200)
}

# Add metadata
for(theSample in names(ctrl.samples.raw)){
  cat("Processing:", theSample, "\n")
  ctrl.samples.seurat[[theSample]] <- AddMetaData(object   = ctrl.samples.seurat[[theSample]], 
                                                  metadata = theSample, 
                                                  col.name = "Patient")
  ctrl.samples.seurat[[theSample]] <- AddMetaData(object   = ctrl.samples.seurat[[theSample]], 
                                                  metadata = "PBMC", 
                                                  col.name = "Tissue")
} 

# Normalise and reduce dimensions
s.genes   <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
for(theSample in names(ctrl.samples.raw)){
  cat("Processing:", theSample, "\n")
  # store mitochondrial percentage in object meta data
  ctrl.samples.seurat[[theSample]] <- PercentageFeatureSet(ctrl.samples.seurat[[theSample]], pattern = "^MT-", col.name = "percent.mt")
  
  # Cell cycle scoring
  ctrl.samples.seurat[[theSample]] <- CellCycleScoring(ctrl.samples.seurat[[theSample]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  # Run sctransform
  ctrl.samples.seurat[[theSample]] <- SCTransform(ctrl.samples.seurat[[theSample]], verbose = T, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  
  # Run PCA and UMAP embedding
  ctrl.samples.seurat[[theSample]] <- RunPCA(ctrl.samples.seurat[[theSample]], verbose = T)
  ctrl.samples.seurat[[theSample]] <- RunUMAP(ctrl.samples.seurat[[theSample]], dims = 1:30, verbose = T)
}

## Integrate samples
# Select features and prep them
features            <- SelectIntegrationFeatures(object.list = ctrl.samples.seurat, nfeatures = 3000)
ctrl.samples.seurat <- PrepSCTIntegration(object.list = ctrl.samples.seurat, anchor.features = features)
ctrl.samples.seurat <- lapply(X = ctrl.samples.seurat, FUN = RunPCA, features = features)

# Find anchors
immune.anchors <- FindIntegrationAnchors(object.list          = ctrl.samples.seurat, 
                                         anchor.features      = features,
                                         normalization.method = "SCT",
                                         dims                 = 1:30, 
                                         reduction            = "rpca",
                                         k.anchor             = 20)

# Integrate the data
ctrl.combined <- IntegrateData(anchorset            = immune.anchors, 
                               normalization.method = "SCT", 
                               dims                 = 1:30)

## Cluster and plot as QC
DefaultAssay(ctrl.combined) <- "integrated"

ctrl.combined <- RunPCA(ctrl.combined, npcs = 30)
ctrl.combined <- RunUMAP(ctrl.combined, reduction = "pca", dims = 1:30)

set.seed(1)
ctrl.combined <- FindNeighbors(ctrl.combined, dims = 1:30, verbose = T, force.recalc = T)
ctrl.combined <- FindClusters(ctrl.combined, verbose = T, resolution = 1, random.seed = 666)

p1 <- DimPlot(ctrl.combined, label = T, pt.size = 3, label.size = 20, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p2 <- DimPlot(ctrl.combined, label = T, pt.size = 3, label.size = 20, shuffle = T, group.by = "Patient") + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
p1 + p2
ggsave("compare_monocytes_to_PBMC_ref/CTRL UMAP clusters.pdf", width = 15, height = 15)

DimPlot(ctrl.combined, label = TRUE, pt.size = 3, label.size = 10, shuffle = T, split.by = "Patient") + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/CTRL UMAP patient split.pdf", width = 15, height = 15)

customFeature(object = ctrl.combined, features = c("CD3E", "CD14", "FCGR3A", "CD79A"), order = T, name = "compare_monocytes_to_PBMC_ref/CTRL feature QC.pdf", width = 15, height = 15)

## Save the control seurat object
saveRDS(ctrl.combined, "Seurat_Objects/control.combined.Wang.seurat.RDS")
#ctrl.combined <- readRDS(file = "Seurat_Objects/control.combined.Wang.seurat.RDS")

#=================================================================================================================================
## Integrate ctrl and patient samples
##================================================================================================================================
# make the list
ctrl_patient.list <- list(ctrl.combined, mye.refined.pbmc_plaque.seurat)
ctrl_patient.list[[1]] <- SCTransform(ctrl_patient.list[[1]], verbose = T, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
DefaultAssay(ctrl_patient.list[[1]]) <- "SCT"
DefaultAssay(ctrl_patient.list[[2]]) <- "SCT"

# Select features and prep them
features          <- SelectIntegrationFeatures(object.list = ctrl_patient.list, nfeatures = 3000)
ctrl_patient.list <- PrepSCTIntegration(object.list = ctrl_patient.list, anchor.features = features)
ctrl_patient.list <- lapply(X = ctrl_patient.list, FUN = RunPCA, features = features)

# Find anchors
immune.anchors <- FindIntegrationAnchors(object.list          = ctrl_patient.list, 
                                         anchor.features      = features,
                                         normalization.method = "SCT",
                                         dims                 = 1:30, 
                                         reduction            = "rpca",
                                         k.anchor             = 20)

ctrl_patient.combined <- IntegrateData(anchorset = immune.anchors,
                               normalization.method = "SCT", 
                               dims                 = 1:30)
DefaultAssay(ctrl_patient.combined) <- "integrated"

# Add metadata to identify which cell came from which experiment
ctrl_patient.combined <- AddMetaData(ctrl_patient.combined, metadata = substr(colnames(ctrl_patient.combined), start = 1, stop = 2), col.name = "source.individual")
source <- ctrl_patient.combined@meta.data$Patient
source[grep("P", source, invert = T)] <- "Control"
source[grep("P", source)]             <- "Patient"
ctrl_patient.combined <- AddMetaData(ctrl_patient.combined, metadata = source, col.name = "source")

#=========================================================================================================================
## Scale, reduce dimensionality, and cluster the integrated data
##========================================================================================================================
# Scale and run embeddings
ctrl_patient.combined <- ScaleData(ctrl_patient.combined)
ctrl_patient.combined <- RunPCA(ctrl_patient.combined, npcs = 30)
ctrl_patient.combined <- RunUMAP(ctrl_patient.combined, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
ctrl_patient.combined <- FindNeighbors(ctrl_patient.combined, reduction = "pca", dims = 1:30, force.recalc =T)
ctrl_patient.combined <- FindClusters( ctrl_patient.combined, resolution = 0.25, random.seed = 666)

# Plot clusters
DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 3, label.size = 10, shuffle = T) + NoLegend() +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP clusters.pdf", width = 15, height = 15)

# PLot origin
DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, group.by = "source")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP source.pdf", width = 15, height = 15)

# Plot split origin
DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, split.by = "source")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP source split.pdf", width = 15, height = 15)

## Save the control_patient seurat object
saveRDS(ctrl_patient.combined, "Seurat_Objects/control_patient.combined.Wang.seurat.RDS")
#ctrl_patient.combined <- readRDS(file = "Seurat_Objects/control_patient.combined.Wang.seurat.RDS")


#=========================================================================================================================
## Map back our refined idendtities and L2 idents from the ref
##========================================================================================================================
## First map to the annotqated PBMC reference set
# Find anchors
anchors <- FindTransferAnchors(reference            = reference,
                               query                = ctrl_patient.combined,
                               normalization.method = "SCT",
                               reference.reduction  = "spca",
                               dims = 1:50
)

# Map to reference
ctrl_patient.combined <- MapQuery(anchorset         = anchors,
                               query                = ctrl_patient.combined,
                               reference            = reference,
                               refdata              = list(celltype.l1   = "celltype.l1",
                                                           celltype.l2   = "celltype.l2",
                                                           predicted_ADT = "ADT"),
                               reference.reduction  = "spca", 
                               reduction.model      = "wnn.umap"
)

# Plot clusters
DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, label.size = 0, shuffle = T, group.by = "predicted.celltype.l1")+
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP clusters ref names L1.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, label.size = 0, shuffle = T, group.by = "predicted.celltype.l2")+
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP clusters ref names L2.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, label.size = 0, shuffle = T, group.by = "source")+
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP clusters source.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, split.by = "source", group.by = "predicted.celltype.l1")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP source split names.pdf", width = 15, height = 15)

## Focus on the myeloid clusters
safetysave.ctrl_patient.combined <- ctrl_patient.combined

# Subset the mono labeled cells
ctrl_patient.combined <- subset(ctrl_patient.combined, predicted.celltype.l1 == "Mono")

# recluster
ctrl_patient.combined <- RunPCA(ctrl_patient.combined, npcs = 30)
ctrl_patient.combined <- RunUMAP(ctrl_patient.combined, reduction = "pca", dims = 1:30)

set.seed(1)
ctrl_patient.combined <- FindNeighbors(ctrl_patient.combined, dims = 1:30, verbose = T, force.recalc = T)
ctrl_patient.combined <- FindClusters(ctrl_patient.combined, verbose = T, resolution = 1, random.seed = 666)

# Plot clusters
DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, label.size = 0, shuffle = T)+
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma clusters.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, label.size = 0, shuffle = T, group.by = "predicted.celltype.l1")+
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma clusters ref names L1.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, label.size = 0, shuffle = T, group.by = "predicted.celltype.l2")+
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma clusters ref names L2.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, label.size = 0, shuffle = T, group.by = "source")+
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma clusters source.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, split.by = "source", group.by = "predicted.celltype.l1")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma source split names.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, split.by = "source")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma source split.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, split.by = "Tissue")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma tissue split.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, split.by = "Tissue", group.by = "source")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma tissue source split.pdf", width = 15, height = 15)


##==================================================================================
##==================================================================================
# Map back sub clustering identities from the patient seurat object
# Pops we want
Mye.idents        <- as.vector(Idents(mye.refined.pbmc_plaque.seurat))
names(Mye.idents) <- names(Idents(mye.refined.pbmc_plaque.seurat))

# Add to the rest of the pops
sub.idents                    <- as.vector(ctrl_patient.combined$predicted.celltype.l2)
names(sub.idents)             <- row.names(ctrl_patient.combined@meta.data)
sub.idents[names(Mye.idents)] <- Mye.idents

# And add to the seurat object
Idents(ctrl_patient.combined) <- sub.idents

# And plot
DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0)  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma clusters proper names.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, split.by = "source")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma source split proper names.pdf", width = 15, height = 15)

DimPlot(ctrl_patient.combined, label = TRUE, pt.size = 2, shuffle = T, label.size = 0, split.by = "Tissue")  +
  theme_pubr(base_size = 14, x.text.angle = 45) +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )
ggsave("compare_monocytes_to_PBMC_ref/UMAP moma tissue split proper names.pdf", width = 15, height = 15)


