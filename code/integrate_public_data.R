#=====================================================================================
## Add in Dib et al. and Gianarelli et al. Macrophages for better resolution and power
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
# Seurat
workdir <- "raw_data/Dib_from_myel_clean_cluster.0.8.dir/"
dib.seurat <- readRDS(paste0(workdir,"begin.rds"))
dib.seurat

# UMAP coords
umap <- read.delim(paste0(workdir,"umap.0.5.tsv.gz"))
head(umap)

# Idents
cl <- read.table(paste0(workdir,"cluster_assignments.tsv.gz"), sep="\t", header=TRUE)
head(cl)
table(cl$cluster_id)

# UMAP coords embedding test
meta <- dib.seurat@meta.data
meta$barcode <- row.names(meta)
meta <- merge(meta, umap, by = 'barcode')
meta <- merge(meta, cl, by = 'barcode')
meta <- meta %>% mutate(cluster_id = factor(cluster_id))
head(meta)

ggplot(meta, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(color=cluster_id), size=1.1) +
  theme_classic() +
  labs(x="UMAP 1", y="UMAP 2") +
  labs(color = "Cluster id")


# Add idents to Seurat object
cl.vec             <- cl$cluster_id
names(cl.vec)      <- cl$barcode
Idents(dib.seurat) <- cl.vec
head(Idents(dib.seurat))

## Prep for integration
dib.seurat <- UpdateSeuratObject(dib.seurat)
dib.seurat <- SCTransform(dib.seurat)
dib.seurat <- RunPCA(dib.seurat)
umap.mat <- as.matrix(meta[,c("UMAP_1", "UMAP_2")])
row.names(umap.mat) <- meta$barcode
dib.seurat[["umap"]] <- CreateDimReducObject(embeddings = umap.mat, assay = "RNA", key = "UMAP_")

# Plot the given populations
customUMAP(object = dib.seurat, pt.size = 2, label = F, title = "Dib et al. data", shuffle = T, seed = 666, plot.width = 10, file.name = "more_macs/dib.UMAP.pdf")
