#======================================================================
## Analyze Interferon signatures in the mac populations
#======================================================================
dir.create("infmac", showWarnings = F)


##==================================================
##==================================================
## Construct signatures using GOBP data 
##==================================================
##==================================================

## Load and preprocess GO genesets
GO.files <- list.files("raw_data/", pattern = "^GO.+INTERFERON.+grp")
geneSets <- list()
for(theGO in GO.files){
  geneSets[[gsub(pattern = ".v2023.2.Hs.grp", replacement = "", x = theGO, fixed = T)]] <- scan(file = paste("raw_data/", theGO, sep = ""), what = "list")
  geneSets[[gsub(pattern = ".v2023.2.Hs.grp", replacement = "", x = theGO, fixed = T)]] <- geneSets[[gsub(pattern = ".v2023.2.Hs.grp", replacement = "", x = theGO, fixed = T)]][4:length(geneSets[[gsub(pattern = ".v2023.2.Hs.grp", replacement = "", x = theGO, fixed = T)]])]
}

##==================================================
## Check for enrichment in the single cell data
# Retrieve the expression matrix from the seurat object
exprMat <- GetAssayData(final.pop.call.integrated.mye.seurat, "counts")
DimPlot(final.pop.call.integrated.mye.seurat) + theme(legend.position = "none")

# Calculate enrichment scores
cells_AUC <- AUCell_run(exprMat, geneSets, BiocParallel::MulticoreParam(5))

# Optimize tresholds
set.seed(333)
par(mfrow=c(2,1)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

cells_assignment$GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON$aucThr$selected
cells_assignment$GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY$aucThr$selected
cells_assignment$GOBP_RESPONSE_TO_TYPE_I_INTERFERON$aucThr$selected
cells_assignment$GOBP_RESPONSE_TO_TYPE_III_INTERFERON$aucThr$selected
cells_assignment$GOBP_TYPE_I_INTERFERON_PRODUCTION$aucThr$selected

length(cells_assignment$GOBP_RESPONSE_TO_TYPE_III_INTERFERON$assignment)


# Add to seurat metadata
inf.sig <- data.frame(row.names = colnames(final.pop.call.integrated.mye.seurat), "inf.sig" = rep("none", length(colnames(final.pop.call.integrated.mye.seurat))))
inf.sig[cells_assignment$GOBP_RESPONSE_TO_TYPE_III_INTERFERON$assignment,  "inf.sig"] <- "IFN-I"
inf.sig[cells_assignment$GOBP_RESPONSE_TO_TYPE_III_INTERFERON$assignment,  "inf.sig"] <- "IFN-II"
inf.sig[cells_assignment$GOBP_RESPONSE_TO_TYPE_III_INTERFERON$assignment,  "inf.sig"] <- "IFN-III"

final.pop.call.integrated.mye.seurat <- AddMetaData(final.pop.call.integrated.mye.seurat, metadata = inf.sig, col.name = "inf.sig")
DimPlot(final.pop.call.integrated.mye.seurat, group.by = "inf.sig")
customUMAP(object = final.pop.call.integrated.mye.seurat, group.by = "inf.sig", pt.size = 2, shuffle = T, legend.pos = "top", seed = 666, file.name = "infmac/GO assignment umap.pdf", cols = c("IFN-III" = "goldenrod", "none" = "grey"))

# Add continuous values
final.pop.call.integrated.mye.seurat <- AddMetaData(object = final.pop.call.integrated.mye.seurat, metadata = as.data.frame(cells_AUC@assays@data$AUC["GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON",]), col.name = "GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON")
final.pop.call.integrated.mye.seurat <- AddMetaData(object = final.pop.call.integrated.mye.seurat, metadata = as.data.frame(cells_AUC@assays@data$AUC["GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY",]), col.name = "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY")
final.pop.call.integrated.mye.seurat <- AddMetaData(object = final.pop.call.integrated.mye.seurat, metadata = as.data.frame(cells_AUC@assays@data$AUC["GOBP_RESPONSE_TO_TYPE_I_INTERFERON",]), col.name = "GOBP_RESPONSE_TO_TYPE_I_INTERFERON")
final.pop.call.integrated.mye.seurat <- AddMetaData(object = final.pop.call.integrated.mye.seurat, metadata = as.data.frame(cells_AUC@assays@data$AUC["GOBP_RESPONSE_TO_TYPE_III_INTERFERON",]), col.name = "GOBP_RESPONSE_TO_TYPE_III_INTERFERON")
final.pop.call.integrated.mye.seurat <- AddMetaData(object = final.pop.call.integrated.mye.seurat, metadata = as.data.frame(cells_AUC@assays@data$AUC["GOBP_TYPE_I_INTERFERON_PRODUCTION",]), col.name = "GOBP_TYPE_I_INTERFERON_PRODUCTION")

# And plot
FeaturePlot(final.pop.call.integrated.mye.seurat, features = "GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON", pt.size = 2)
FeaturePlot(final.pop.call.integrated.mye.seurat, features = "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY", pt.size = 2)
FeaturePlot(final.pop.call.integrated.mye.seurat, features = "GOBP_RESPONSE_TO_TYPE_I_INTERFERON", pt.size = 2)
FeaturePlot(final.pop.call.integrated.mye.seurat, features = "GOBP_RESPONSE_TO_TYPE_III_INTERFERON", pt.size = 2)
FeaturePlot(final.pop.call.integrated.mye.seurat, features = "GOBP_TYPE_I_INTERFERON_PRODUCTION", pt.size = 2)

customFeature(object = final.pop.call.integrated.mye.seurat, features = "GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON", pt.size = 4, name = "infmac/GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON score umap.pdf", cols = c("grey", "red3"))
customFeature(object = final.pop.call.integrated.mye.seurat, features = "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY", pt.size = 4, name = "infmac/GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY score umap.pdf", cols = c("grey", "red3"))
customFeature(object = final.pop.call.integrated.mye.seurat, features = "GOBP_RESPONSE_TO_TYPE_I_INTERFERON", pt.size = 4, name = "infmac/GOBP_RESPONSE_TO_TYPE_I_INTERFERON score umap.pdf", cols = c("grey", "red3"))
customFeature(object = final.pop.call.integrated.mye.seurat, features = "GOBP_RESPONSE_TO_TYPE_III_INTERFERON", pt.size = 4, name = "infmac/GOBP_RESPONSE_TO_TYPE_III_INTERFERON score umap.pdf", cols = c("grey", "red3"))
customFeature(object = final.pop.call.integrated.mye.seurat, features = "GOBP_TYPE_I_INTERFERON_PRODUCTION", pt.size = 4, name = "infmac/GOBP_TYPE_I_INTERFERON_PRODUCTION score umap.pdf", cols = c("grey", "red3"))


##==================================================
##==================================================
## Construct signatures using custom IFNCI signature
##==================================================
##==================================================
## Retrieve IFNIC signature from the IFN high population
# Get Interferon cluster marker genes
IFN.pop     <- as.vector(unique(Idents(final.pop.call.integrated.mye.seurat))[grep(pattern = "Interferon", x = unique(Idents(final.pop.call.integrated.mye.seurat)))])
IFN.markers <- FindMarkers(final.pop.call.integrated.mye.seurat, ident.1 = IFN.pop, assay = "RNA", min.pct = 0.25, min.diff.pct = 0.25, only.pos = T)
IFN.markers <- row.names(IFN.markers[IFN.markers$p_val_adj < 0.05,])

# Compare with the GO sets
all.IFN.GO <- unlist(geneSets)
all.IFN.GO <- unique(all.IFN.GO)

IFN.markers.GO <- IFN.markers[IFN.markers %in% all.IFN.GO]

geneSets <- list("IFNIC.full" = IFN.markers, "IFNIC.GO" = IFN.markers.GO)

##==================================================
## Check for enrichment in the single cell data
# Retrieve the expression matrix from the seurat object
exprMat <- GetAssayData(final.pop.call.integrated.mye.seurat, "counts")
DimPlot(final.pop.call.integrated.mye.seurat) + theme(legend.position = "none")

# Calculate enrichment scores
cells_AUC <- AUCell_run(exprMat, geneSets, BiocParallel::MulticoreParam(5))

# Optimize tresholds
set.seed(333)
par(mfrow=c(2,1)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

cells_assignment$IFNIC.full$aucThr$selected
cells_assignment$IFNIC.GO$aucThr$selected
length(cells_assignment$IFNIC.full$assignment)
length(cells_assignment$IFNIC.GO$assignment)


# Add discrete assignments to seurat metadata
inf.sig <- data.frame(row.names = colnames(final.pop.call.integrated.mye.seurat), "inf.sig" = rep("none", length(colnames(final.pop.call.integrated.mye.seurat))))
inf.sig[cells_assignment$IFNIC.full$assignment,  "inf.sig"] <- "IFNIC.full"
inf.sig[cells_assignment$IFNIC.GO$assignment,  "inf.sig"] <- "IFNIC.GO"
final.pop.call.integrated.mye.seurat <- AddMetaData(final.pop.call.integrated.mye.seurat, metadata = inf.sig, col.name = "inf.sig")
                                                     
# And plot
DimPlot(final.pop.call.integrated.mye.seurat, group.by = "inf.sig")
customUMAP(object = final.pop.call.integrated.mye.seurat, group.by = "inf.sig", pt.size = 2, shuffle = T, legend.pos = "top", seed = 666, file.name = "infmac/IFNIC assignment umap.pdf", cols = c("IFNIC.full" = "goldenrod", "IFNIC.GO" = "orangered", "none" = "grey"))
customUMAP(object = final.pop.call.integrated.mye.seurat, pt.size = 2, shuffle = T, legend.pos = "none", seed = 666, file.name = "infmac/Ref umap.pdf", cols = M.int_refined.pop.colors)
customUMAP(object = final.pop.call.integrated.mye.seurat, pt.size = 2, shuffle = T, legend.pos = "right", seed = 666, file.name = "infmac/Ref umap legend.pdf", cols = M.int_refined.pop.colors)


# Add continuous values
final.pop.call.integrated.mye.seurat <- AddMetaData(object = final.pop.call.integrated.mye.seurat, metadata = as.data.frame(cells_AUC@assays@data$AUC["IFNIC.full",]), col.name = "IFNIC.full")
final.pop.call.integrated.mye.seurat <- AddMetaData(object = final.pop.call.integrated.mye.seurat, metadata = as.data.frame(cells_AUC@assays@data$AUC["IFNIC.GO",]), col.name = "IFNIC.GO")

# And plot
FeaturePlot(final.pop.call.integrated.mye.seurat, features = "IFNIC.full", pt.size = 2)
FeaturePlot(final.pop.call.integrated.mye.seurat, features = "IFNIC.GO", pt.size = 2)
customFeature(object = final.pop.call.integrated.mye.seurat, features = "IFNIC.full", pt.size = 4, name = "infmac/IFNIC.full score umap.pdf", cols = c("grey", "red3"))
customFeature(object = final.pop.call.integrated.mye.seurat, features = "IFNIC.GO", pt.size = 4, name = "infmac/IFNIC.GO score umap.pdf", cols = c("grey", "red3"))


