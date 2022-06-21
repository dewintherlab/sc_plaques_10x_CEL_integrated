#=========================================================================================================================
## Integrate  all cells with the 43 patient CEL-seq myeloid cells
##========================================================================================================================
dir.create("full_43p_celseq_integration/First round of clustering/clusters", showWarnings = F, recursive = T)

## Load the data
full.43p.seurat <- readRDS(file = "raw_data/main.seurat.complete.RDS")
full.43p.seurat
unique(Idents(full.43p.seurat))


## Prep for integration
# Select the corect assay
DefaultAssay(full.43p.seurat) <- "RNA"
DefaultAssay(patients.seurat) <- "RNA"

# Make a list of our objects to integrate
integrated.full.seurat <- list(patients.seurat, full.43p.seurat)

# Define var features
integrated.full.seurat <- lapply(X = integrated.full.seurat, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = integrated.full.seurat)

# Scale data and run PCA using the selected features
integrated.full.seurat <- lapply(X = integrated.full.seurat, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

## Integrate the data
## Use RPCA  to speed things up.
# First select anchors for integration (by RPCA reduction)
anchors <- FindIntegrationAnchors(object.list = integrated.full.seurat, anchor.features = features, reduction = "rpca")

# And do the actual integrating
integrated.full.seurat <- IntegrateData(anchorset = anchors)

# Show stats
integrated.full.seurat

#=========================================================================================================================
## Update metadata
head(integrated.full.seurat@meta.data)
tail(integrated.full.seurat@meta.data)

# Fix tissue column from 43p to 'plaque'
integrated.full.seurat@meta.data[is.na(integrated.full.seurat$Tissue),"Tissue"] <- "plaque"

# Fix sex column from 10X to 'male'
integrated.full.seurat@meta.data[is.na(integrated.full.seurat$Sex),"Sex"] <- "male"

# Add a new column 'Method' to distinguish the two datasets
method <- rep("CEL-seq", dim(integrated.full.seurat)[2])
method[grep("^P", integrated.full.seurat$Patient)] <- "10X"
integrated.full.seurat <- AddMetaData(integrated.full.seurat, method, col.name = "Method")

## Add a new column 'original.refined.ident' to capture the original idents of the unintegrated objects.
# Set up destination dataframe
ref.ident  <- data.frame(row.names = row.names(integrated.full.seurat@meta.data))

# Fetch idents and stick 'em in a dataframe for merging
idents.43p <- data.frame(refined.ident = Idents(full.43p.seurat))
idents.10X <- data.frame(refined.ident = Idents(patients.seurat))
idents.com <- rbind(idents.43p, idents.10X)

# Merge the frames together
ref.ident <- merge(ref.ident, idents.com, by = 0)
row.names(ref.ident) <- ref.ident$Row.names
ref.ident$Row.names <- NULL
tail(ref.ident)

# And add the metadata
integrated.full.seurat <- AddMetaData(integrated.full.seurat, metadata = ref.ident, col.name = "original.refined.ident")
                                      

#=========================================================================================================================
## Run the embeddings, cluster, and have a look
##========================================================================================================================
# Run dimensionality reductions
integrated.full.seurat <- ScaleData(object = integrated.full.seurat, verbose = F)
integrated.full.seurat <- RunPCA(   object = integrated.full.seurat, verbose = F)
integrated.full.seurat <- RunUMAP(  object = integrated.full.seurat, reduction = "pca", dims = 1:30)

# Cluster
set.seed(1)
integrated.full.seurat <- FindNeighbors(object = integrated.full.seurat, reduction = "pca", dims = 1:30, force.recalc = T, verbose = T)
integrated.full.seurat <- FindClusters( object = integrated.full.seurat, resolution = 1.25, verbose = T, random.seed = 666)

## Visualize
# Plot Clusters
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", 
           pt.size    = 2,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration intitial clustering UMAP.pdf", 
           plot.width = 15)

# Plot Method
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Method",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration Method UMAP.pdf", 
           plot.width = 15)

# Plot Patients
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Patient",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration Patient UMAP.pdf", 
           plot.width = 15)

# Plot Tissue
customUMAP(object     = integrated.full.seurat, shuffle = T,
           title      = "Initial clustering", 
           legend.pos = "right", group.by = "Tissue",
           pt.size    = 1,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration Tissue UMAP.pdf", 
           plot.width = 15)



## Plot the original 10X and 43p idents
customUMAP(object     = integrated.full.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "Source clusters", legend.pos = "none", 
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration original cluster UMAP.pdf", 
           plot.width = 15)

p1 <- DimPlot(integrated.full.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, pt.size = 3) +
  ggtitle("Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.full.seurat, reduction = "umap", label = TRUE, repel = TRUE, shuffle = T, pt.size = 3) +
  ggtitle("Integrated clusters") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("full_43p_celseq_integration/First round of clustering/Integration original vs. integrated cluster UMAP.pdf", width = 20, height = 10)

# Highlight each source cluster
dir.create("full_43p_celseq_integration/First round of clustering/Original source cluster highlights", showWarnings = F, recursive = T)
for(theIdent in unique(integrated.full.seurat$original.refined.ident)){
  customUMAP(object          = integrated.full.seurat, reduction = "umap", label = F, shuffle = T, title = theIdent, legend.pos = "none", 
             file.name       = paste("full_43p_celseq_integration/First round of clustering/Original source cluster highlights/", theIdent, " UMAP.pdf", sep = ""),
             sizes.highlight = 4,
             cells.highlight = row.names(ref.ident)[which(ref.ident$refined.ident == theIdent)])
}

## Color only one source at a time
# 43p colored
cols.43p        <- colorSpacer(startcolor = "red", endcolor = "blue", steps = length(unique(Idents(subset(integrated.full.seurat, Method != "10X")))))
names(cols.43p) <- unique(Idents(subset(integrated.full.seurat, Method != "10X")))
customUMAP(object     = integrated.full.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "CEL-seq Source population identities", legend.pos = "none", 
           pt.size    = 5,
           cols       = cols.43p,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration original 43p cluster UMAP.pdf", 
           plot.width = 15)

# 10X colored
cols.10X        <- colorSpacer(startcolor = "yellow", endcolor = "green", steps = length(unique(Idents(subset(integrated.full.seurat, Method == "10X")))))
names(cols.10X) <- unique(Idents(subset(integrated.full.seurat, Method == "10X")))
customUMAP(object     = integrated.full.seurat, reduction = "umap", group.by = "original.refined.ident", shuffle = T, title = "10X Source population identities", legend.pos = "none", 
           pt.size    = 5,
           cols       = cols.10X,
           file.name  = "full_43p_celseq_integration/First round of clustering/Integration original 10X cluster UMAP.pdf", 
           plot.width = 15)

# Double plot
p1 <- DimPlot(integrated.full.seurat, reduction = "umap", cols = cols.43p, group.by = "original.refined.ident", shuffle = T, pt.size = 5) +
  ggtitle("CEL-seq Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p2 <- DimPlot(integrated.full.seurat, reduction = "umap", cols = cols.10X, group.by = "original.refined.ident", shuffle = T, pt.size = 5) +
  ggtitle("10X Source population identities") +
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "none") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold")
  )

p1 + p2
ggsave("full_43p_celseq_integration/First round of clustering/Integration 43p vs. 10X clusters UMAP.pdf", width = 20, height = 10)

## Save the integrated seurat object
saveRDS(integrated.full.seurat, "Seurat_Objects/full.43p_10X.integrated.seurat.RDS")
#integrated.full.seurat <- readRDS(file = "Seurat_Objects/full.43p_10X.integrated.seurat.RDS")

## Incorporate our refined monocyte and macrophage identities
# Get myeloid idents 
mye.idents <- c(Idents(integrated.mye.seurat), Idents(pbmc.monocytes.seurat))
mye.cells  <- names(mye.idents)

# Set myeloid idents
Idents(integrated.full.seurat, cells = mye.cells) <- mye.idents

# Unset all others, for now
Idents(integrated.full.seurat, cells = colnames(integrated.full.seurat)[!colnames(integrated.full.seurat) %in% mye.cells]) <- "Other"

# Update colors
tmp.col                 <- colorSpacer(startcolor = "grey", steps = 1)
names(tmp.col)          <- "Other"
full.refined.pop.colors <- c(M.int_refined.pop.colors, tmp.col)

customUMAP(object = integrated.full.seurat, legend.pos = "right", plot.width = 20, cols = full.refined.pop.colors,
           file.name = "full_43p_celseq_integration/First round of clustering/Refined mye pops UMAP.pdf")


