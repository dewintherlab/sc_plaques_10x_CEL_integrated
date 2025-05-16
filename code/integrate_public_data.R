#=====================================================================================
## Add in Dib et al. and Horstmann et al. Macrophages for better resolution and power
#=====================================================================================
## Create dir
dir.create("more_macs", showWarnings = F)

##================##
## Load our data. ##
##================##
final.pop.call.integrated.mac.seurat <- readRDS(file = "Seurat_Objects/final.pop.call.integrated.mac.seurat.RDS")
archetype.colors                     <- readRDS(file = "Seurat_Objects/archetype.colors.RDS")
M.int_refined.pop.colors             <- readRDS(file = "Seurat_Objects/M.int_refined.pop.colors.RDS")

# Plot a UMAP for reference
customUMAP(object = final.pop.call.integrated.mac.seurat, pt.size = 3, label = F, title = "Prange et al. data", shuffle = T, seed = 666, plot.width = 10, file.name = "more_macs/prange.UMAP.pdf")

##===========================##
## Load the Dib et al. data. ##
##===========================##
## Load the objects
## COMMENTED OUT SINCE AFTER THE FIRST TIME WE SAVED THE OBJECT AND IT'S QUICKER TO JUST READ THAT IN


# Seurat
# workdir <- "raw_data/Dib_from_myel_clean_cluster.0.8.dir/"
# dib.seurat <- readRDS(paste0(workdir,"begin.rds"))
# dib.seurat
# 
# # UMAP coords
# umap <- read.delim(paste0(workdir,"umap.0.5.tsv.gz"))
# head(umap)
# 
# # Idents
# cl <- read.table(paste0(workdir,"cluster_assignments.tsv.gz"), sep="\t", header=TRUE)
# head(cl)
# table(cl$cluster_id)
# 
# # UMAP coords embedding test
# meta <- dib.seurat@meta.data
# meta$barcode <- row.names(meta)
# meta <- merge(meta, umap, by = 'barcode')
# meta <- merge(meta, cl, by = 'barcode')
# meta <- meta %>% mutate(cluster_id = factor(cluster_id))
# head(meta)
# 
# ggplot(meta, aes(x=UMAP_1, y=UMAP_2)) +
#   geom_point(aes(color=cluster_id), size=1.1) +
#   theme_classic() +
#   labs(x="UMAP 1", y="UMAP 2") +
#   labs(color = "Cluster id")
# 
# 
# # Add idents to Seurat object
# cl.vec             <- cl$cluster_id
# names(cl.vec)      <- cl$barcode
# Idents(dib.seurat) <- cl.vec
# head(Idents(dib.seurat))
# 
# ## Prep for integration
# dib.seurat <- UpdateSeuratObject(dib.seurat)
# dib.seurat <- SCTransform(dib.seurat)
# dib.seurat <- RunPCA(dib.seurat)
# umap.mat <- as.matrix(meta[,c("UMAP_1", "UMAP_2")])
# row.names(umap.mat) <- meta$barcode
# dib.seurat[["umap"]] <- CreateDimReducObject(embeddings = umap.mat, assay = "RNA", key = "UMAP_")
# 
# saveRDS(dib.seurat, "Seurat_Objects/dib.seurat.RDS")

dib.seurat <- readRDS(file = "Seurat_Objects/dib.seurat.RDS")

# Plot the given populations
customUMAP(object = dib.seurat, pt.size = 2, label = F, title = "Dib et al. data", shuffle = T, seed = 666, plot.width = 10, file.name = "more_macs/dib.UMAP.pdf")

##=================================##
## Load the Horstmann et al. data. ##
##=================================##
## Load the samples
# Fetch sample dirs
all.sample.dirs <- list.files(path = "raw_data/Horstmann", include.dirs = T, pattern = "P", full.names = T)

# Load raw data
samples.raw <- list()
for(theSample in all.sample.dirs){
  cat("Loading sample:", basename(theSample), "...\n")
  
  # Load from disk
  samples.raw[[basename(theSample)]] <- Read10X(data.dir = theSample)
}

# Turn into seurat objects
samples.seurat <- list()
for(theSample in names(samples.raw)){
  cat("Processing:", theSample, "\n")
  samples.seurat[[theSample]] <- CreateSeuratObject(counts       = samples.raw[[theSample]],
                                                    project      = "scTCRseq", 
                                                    min.cells    = 3,
                                                    min.features = 200)
}

# Clean up
rm(samples.raw)
rm(all.sample.dirs)

## Normalize and cluster
# run sctransform
samples.seurat <- lapply(samples.seurat, function(x)SCTransform(x, verbose = T))

# Run PCA and UMAP embedding
samples.seurat <- lapply(samples.seurat, function(x)RunPCA( x, verbose = T))
samples.seurat <- lapply(samples.seurat, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Cluster
set.seed(1)
samples.seurat <- lapply(samples.seurat, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
samples.seurat <- lapply(samples.seurat, function(x)FindClusters( x, verbose = T, resolution = 1, random.seed = 666))

# QC
samples.seurat <- lapply(samples.seurat, function(x)PercentageFeatureSet(x, pattern = "^MT-",     col.name = "percent.mt"))
samples.seurat <- lapply(samples.seurat, function(x)PercentageFeatureSet(x, pattern = "KCNQ1OT1", col.name = "percent.KCNQ1OT1"))
samples.seurat <- lapply(samples.seurat, function(x)PercentageFeatureSet(x, pattern = "UGDH-AS1", col.name = "percent.UGDH.AS1"))
samples.seurat <- lapply(samples.seurat, function(x)PercentageFeatureSet(x, pattern = "GHET1",    col.name = "percent.GHET1"))

for(theSample in names(samples.seurat)){
  VlnPlot(object = samples.seurat[[theSample]], features = c("percent.mt", "percent.KCNQ1OT1", "percent.UGDH.AS1", "percent.GHET1"))
  ggsave(paste("more_macs/Horstmann_QC_", theSample, " Percentage mmt and Apoptotic markers.pdf", sep = ""))
}

# Keep only good cells
samples.seurat <- lapply(samples.seurat, function(x)subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 ))

## Rerun normalisation while regressing out mitochondrial reads
vars.to.regress <- c("percent.mt")

# run sctransform
samples.seurat <- lapply(samples.seurat, function(x)SCTransform(x, vars.to.regress = vars.to.regress, verbose = T))

# Run PCA and UMAP embedding
samples.seurat <- lapply(samples.seurat, function(x)RunPCA( x, verbose = T))
samples.seurat <- lapply(samples.seurat, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Cluster
set.seed(1)
samples.seurat <- lapply(samples.seurat, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
samples.seurat <- lapply(samples.seurat, function(x)FindClusters( x, verbose = T, resolution = 1, random.seed = 666))

## Integrate
# Prep for integration
features        <- SelectIntegrationFeatures(object.list         = samples.seurat, nfeatures = 3000)
samples.seurat  <- PrepSCTIntegration(       object.list         = samples.seurat, anchor.features = features)
anchors         <- FindIntegrationAnchors(   object.list         = samples.seurat, 
                                             normalization.method = "SCT",
                                             anchor.features      = features, 
                                             dims                 = 1:30, 
                                             reduction            = "cca", 
                                             k.anchor             = 20)

# Save it
save(anchors, samples.seurat, file = "Seurat_Objects/horstmann.seurat_list.anchors.Rdata")

# Integrate the samples
samples.seurat <- IntegrateData(anchorset             = anchors, 
                                normalization.method = "SCT", 
                                dims                 = 1:30)

# Run PCA and UMAP embedding
samples.seurat <- lapply(samples.seurat, function(x)RunPCA( x, verbose = T))
samples.seurat <- lapply(samples.seurat, function(x)RunUMAP(x, dims = 1:30, verbose = T))

# Cluster
set.seed(1)
samples.seurat <- lapply(samples.seurat, function(x)FindNeighbors(x, dims = 1:30, verbose = T, force.recalc = T))
samples.seurat <- lapply(samples.seurat, function(x)FindClusters( x, verbose = T, resolution = 1, random.seed = 666))

## Scale all genes from RNA assay for diff gene calling and visualisation
# Set default to RNA
DefaultAssay(samples.seurat) <- "RNA"

# Normalize
samples.seurat <- NormalizeData(samples.seurat, normalization.method = "LogNormalize")

# Scale
samples.seurat <- ScaleData(samples.seurat, features = row.names(samples.seurat))

# Look at some genes to determine where the macrophages are
bunchOfCustomPlots(object    = samples.seurat,
                   features  = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                 "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                 "CD79A", "CD86", "CD34", "KIT", "OLR1", 
                                 "FCGR3A", "SELL", "FOLR2", "TREM1", "PLIN2"), 
                   ncol      = 5,
                   name      = "more_macs/Horstmann Cell type defining genes",
                   Vln.width = 25)
