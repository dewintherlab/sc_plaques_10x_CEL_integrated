## ssGSEA analysis with ESCAPE
dir.create("ssGSEA_ssGSEA_results", showWarnings = F)

# Select genesets of interest from msigdb
as.data.frame(msigdbr::msigdbr_collections())
gsets <- getGeneSets(library = c("H", "C2","C5", "C8"))

# Remove CGP subcat
cgp   <- unique(msigdbr(category = "C2", subcategory = "CGP")$gs_name)
gsets <- gsets[!names(gsets) %in% cgp,]

# Keep only GO biological process
notGOBP  <- c(unique(msigdbr(category = "C5", subcategory = c("GO:CC"))$gs_name),
              unique(msigdbr(category = "C5", subcategory = c("GO:MF"))$gs_name),
              unique(msigdbr(category = "C5", subcategory = c("HPO"))$gs_name))
gsets <- gsets[!names(gsets) %in% notGOBP,]

# Alternatively, don't use GO at all (GO, SCHMO)
gsets <- getGeneSets(library = c("H", "C2", "C8"))

# Run escape enrichment
ssGSEA.integrated.mye.seurat <- enrichIt(obj = integrated.mye.seurat, gene.sets = gsets, groups = 1000, cores = 10)
integrated.mye.seurat <- Seurat::AddMetaData(integrated.mye.seurat, ssGSEA.integrated.mye.seurat)

## Save the objects
saveRDS(ssGSEA.integrated.mye.seurat, "Seurat_Objects/ssGSEA_data.43p_10X.integrated.RDS")
#ssGSEA.integrated.mye.seurat <- readRDS(file = "Seurat_Objects/ssGSEA_data.43p_10X.integrated.RDS")

saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#ssGSEA.integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")


## Visualise
# Set up a color palette
hmcolors <- colorRampPalette(c("#0348A6", "#7AC5FF", "#C6FDEC", "#FFB433", "#FF4B20"))

## Heatmap
# Subset highest enrichments
# Base cutoff value on trial and error to get a decent number for a heatmap
cutoff <- 0.49
subset.ssGSEA.sets <- apply(ssGSEA.integrated.mye.seurat,2,function(x) max(abs(x))) > cutoff
subset.ssGSEA.sets <- names(subset.ssGSEA.sets[subset.ssGSEA.sets == 1])

# Check the number
length(subset.ssGSEA.sets)

# Plot it
integrated.mye.seurat <- Seurat::AddMetaData(integrated.mye.seurat, data.frame(cellmethod=Idents(integrated.mye.seurat)))
dittoHeatmap(integrated.mye.seurat, metas = subset.ssGSEA.sets, genes = NULL, annot.by = "cellmethod",
             fontsize = 7,
             cluster_cols = F,
             heatmap.colors = hmcolors(50), filename = "ssGSEA_ssGSEA_results/ssGSEA_top_heatmap.pdf", width=15, height =10)


ES2 <- data.frame(ssGSEA.integrated.mye.seurat, Idents(integrated.mye.seurat), integrated.mye.seurat$Method)
colnames(ES2)[ncol(ES2)-1] <- "cluster"
colnames(ES2)[ncol(ES2)] <- "method"

PCA <- performPCA(enriched = ES2, groups = "cluster")
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = F)

masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)
ggsave("ssGSEA_results/ssGSEA.master_PCA.pdf")

ridgeEnrichment(ES2, gene.set = "HALLMARK_INFLAMMATORY_RESPONSE", facet = "cluster", group = "method", add.rug = TRUE, colors = c("salmon", "dodgerblue"))
ggsave("ssGSEA_results/HALLMARK_INFLAMMATORY_RESPONSE.pdf", width=15, height = 5)

ridgeEnrichment(ES2, gene.set = "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE", facet = "cluster", group = "method", add.rug = TRUE, colors = c("salmon", "dodgerblue"))
ggsave("ssGSEA_results/GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.pdf", width=15, height = 5)

ridgeEnrichment(ES2, gene.set = "GOBP_REGULATION_OF_IMMUNE_RESPONSE", group = "method", facet = "cluster", add.rug = TRUE)
ggsave("ssGSEA_results/GOBP_REGULATION_OF_IMMUNE_RESPONSE.pdf", width=10, height = 5)

ridgeEnrichment(ES2, gene.set = "GOBP_ADAPTIVE_IMMUNE_RESPONSE", group = "method", facet = "cluster", add.rug = TRUE)
ggsave("ssGSEA_results/GOBP_ADAPTIVE_IMMUNE_RESPONSE.pdf", width=10, height = 5)

ridgeEnrichment(ES2, gene.set = "HALLMARK_INFLAMMATORY_RESPONSE", group = "cluster",  add.rug = TRUE)
ggsave("ssGSEA_results/HALLMARK_INFLAMMATORY_RESPONSE.per_cluster.pdf", width=15, height = 10)

ridgeEnrichment(ES2, gene.set = "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE", group = "cluster",  add.rug = TRUE)
ggsave("ssGSEA_results/GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.per_cluster.pdf", width=15, height = 10)

ridgeEnrichment(ES2, gene.set = "GOBP_REGULATION_OF_IMMUNE_RESPONSE", group = "cluster",  add.rug = TRUE)
ggsave("ssGSEA_results/GOBP_REGULATION_OF_IMMUNE_RESPONSE.per_cluster.pdf", width=10, height = 10)

ridgeEnrichment(ES2, gene.set = "GOBP_ADAPTIVE_IMMUNE_RESPONSE", group = "cluster",  add.rug = TRUE)
ggsave("ssGSEA_results/GOBP_ADAPTIVE_IMMUNE_RESPONSE.per_cluster.pdf", width=10, height = 10)


ridgeEnrichment(ES2, gene.set = "GOBP_REGULATION_OF_LEUKOCYTE_MIGRATION", group = "cluster",  add.rug = TRUE, colors = M.int_refined.pop.colors[4:13])
ggsave("ssGSEA_results/GOBP_REGULATION_OF_LEUKOCYTE_MIGRATION_per_cluster.pdf", width=10, height = 10)

ridgeEnrichment(ES2, gene.set = "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION", group = "cluster",  add.rug = TRUE, colors = M.int_refined.pop.colors[4:13])
ggsave("ssGSEA_results/KEGG LEUKOCYTE TRANSENDOTHELIAL MIGRATION_per_cluster.pdf", width=10, height = 10)

ridgeEnrichment(ES2, gene.set = "GOBP_LEUKOCYTE_MIGRATION", group = "cluster",  add.rug = TRUE, colors = M.int_refined.pop.colors[4:13])
ggsave("ssGSEA_results/GOBP_LEUKOCYT_MIGRATION_per_cluster.pdf", width=10, height = 10)

ridgeEnrichment(ES2, gene.set = "GOBP_LEUKOCYTE_CELL_CELL_ADHESION", group = "cluster", add.rug = TRUE, colors = M.int_refined.pop.colors[4:13])
ggsave("ssGSEA_results/GOBP_LEUKOCYTE_CELL_CELL_ADHESION_per_cluster.pdf", width=10, height = 10)


## Plot top 10 pathways per population
for(thePop in names(ont.M)){
  dir.create(paste("ssGSEA_results/", thePop, " top 10 pathways", sep = ""), recursive = T, showWarnings = F)
  cat(paste("Working on ", thePop, "...\n", sep = ""))
  for(thePath in head(ont.M[[thePop]]$Pathways, n = 10)$name){
    cat(paste("\tPlotting: ", thePath, "...\n", sep = ""))
    
    if(!thePath %in% colnames(ES2)){
      cat(paste("\t\tPathway not found!\n", sep = ""))
    }else{
     ridgeEnrichment(ES2, gene.set = thePath, group = "cluster",  add.rug = TRUE, colors = M.int_refined.pop.colors[4:13])
      ggsave(paste("ssGSEA_results/", thePop, " top 10 pathways/", thePath, "_per_cluster.pdf", sep = ""), width=10, height = 10)
    }
  }
}


## Calculate enrichment of marker genes
# Get significant markers
marker.gene.sets <- list()
for(thePop in names(ont.M)){
  marker.gene.sets[[thePop]] <- subset(integrated.mye.seurat.markers, cluster == thePop & p_val_adj < 0.1)$gene
}

# Get enrichment scores
ssGSEA.markers.integrated.mye.seurat <- enrichIt(obj = integrated.mye.seurat, gene.sets = marker.gene.sets, groups = 1000, cores = 10)
integrated.mye.seurat <- Seurat::AddMetaData(integrated.mye.seurat, ssGSEA.markers.integrated.mye.seurat)


## Save the objects
saveRDS(ssGSEA.markers.integrated.mye.seurat, "Seurat_Objects/ssGSEA_markers.data.43p_10X.integrated.RDS")
#ssGSEA.markers.integrated.mye.seurat <- readRDS(file = "Seurat_Objects/ssGSEA_markers.data.43p_10X.integrated.RDS")

saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#ssGSEA.integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")


## Plot it
# Match color names (mangled by forbidden chars in colnames)
ggsea.marker.colors   <- M.int_refined.pop.colors[4:13]
names(ggsea.marker.colors)   <- gsub("+", ".", names(ggsea.marker.colors), fixed = T)
names(ggsea.marker.colors)   <- gsub("-", ".", names(ggsea.marker.colors), fixed = T)
names(ggsea.marker.colors)   <- gsub(" ", ".", names(ggsea.marker.colors), fixed = T)
ggsea.marker.colors <- ggsea.marker.colors[colnames(integrated.mye.seurat@meta.data)[grep("CD14", colnames(integrated.mye.seurat@meta.data))]]

ggsea.marker.idents            <- data.frame(cellmethod=Idents(integrated.mye.seurat))
ggsea.marker.idents$cellmethod <- gsub("+", ".", ggsea.marker.idents$cellmethod, fixed = T)
ggsea.marker.idents$cellmethod <- gsub("-", ".", ggsea.marker.idents$cellmethod, fixed = T)
ggsea.marker.idents$cellmethod <- gsub(" ", ".", ggsea.marker.idents$cellmethod, fixed = T)
integrated.mye.seurat          <- Seurat::AddMetaData(integrated.mye.seurat, ggsea.marker.idents)

unique(integrated.mye.seurat$cellmethod)
names(ggsea.marker.colors)

# Plot goruped by ident
dittoHeatmap(integrated.mye.seurat, metas =  colnames(integrated.mye.seurat@meta.data)[grep("CD14", colnames(integrated.mye.seurat@meta.data))], 
             annot.by = "cellmethod",
             annot.colors = ggsea.marker.colors,
             genes = NULL,
             fontsize = 7,
             cluster_cols = F,
             heatmap.colors = hmcolors(50), filename = "ssGSEA_results/ssGSEA_marker_heatmap.pdf", width=20, height =5)

# Plot cells clustered
dittoHeatmap(integrated.mye.seurat, metas =  colnames(integrated.mye.seurat@meta.data)[grep("CD14", colnames(integrated.mye.seurat@meta.data))], 
             annot.by = "cellmethod",
             annot.colors = ggsea.marker.colors,
             genes = NULL,
             fontsize = 7,
             cluster_cols = T,
             heatmap.colors = hmcolors(50), filename = "ssGSEA_results/ssGSEA_marker_heatmap_cells_clustered.pdf", width=20, height =5)


ES2 <- data.frame(ssGSEA.markers.integrated.mye.seurat, Idents(integrated.mye.seurat), integrated.mye.seurat$Method)
colnames(ES2)[ncol(ES2)-1] <- "cluster"
colnames(ES2)[ncol(ES2)] <- "method"

PCA <- performPCA(enriched = ES2, groups = "cluster")
pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = F)

masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", top.contribution = 10)
ggsave("ssGSEA_results/ssGSEA.marker.master_PCA.pdf")

for(thePath in names(ggsea.marker.colors)){
  cat(paste("Plotting: ", thePath, "...\n", sep = ""))
  if(!thePath %in% colnames(ES2)){
    cat(paste("\t\tPathway not found!\n", sep = ""))
  }else{
    ridgeEnrichment(ES2, gene.set = thePath, group = "cluster",  add.rug = TRUE, colors = M.int_refined.pop.colors[4:13])
    ggsave(paste("ssGSEA_results/marker_", thePath, "_per_cluster.pdf", sep = ""), width=10, height = 10)
  }
}


