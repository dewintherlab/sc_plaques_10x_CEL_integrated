#======================================================================
## Final mac population call and plotting
#======================================================================
dir.create("final_mac_pops/clusters", showWarnings = F, recursive = T)

## First do the mono included object
## Merge the two barely different inf pops (CD16- and S100A8+) into a fresh seurat object to keep things safe
# Check the old idents
unique(Idents(from_full.integrated.mye.seurat))

# And rename
final.pop.call.from_full.integrated.mye.seurat <- RenameIdents(object = from_full.integrated.mye.seurat, "CD14+IL1B+SELL+CD16- Antigen-presenting Inflammatory Monocyte-derived Macrophages" = "CD14+IL1B+SELL+S100A8+ Migrating Inflammatory Monocyte-derived Macrophages")

# Check the new idents
unique(Idents(final.pop.call.from_full.integrated.mye.seurat))

# Update the UMAP
final.pop.call.from_full.integrated.mye.seurat <- RunUMAP(final.pop.call.from_full.integrated.mye.seurat, dims = 1:50, return.model = T)

# Draw the new UMAP
customUMAP(object     = final.pop.call.from_full.integrated.mye.seurat, 
           pt.size    = 3, 
           label      = F, 
           cols       = M.int_refined.pop.colors, 
           shuffle    = T, 
           legend.pos = "right",
           file.name  = "final_mac_pops/population UMAP plus mono.pdf", 
           plot.width = 15)


## Check this one strange mono cluster that popped up on rerunning the embedding
plot         <- DimPlot(final.pop.call.from_full.integrated.mye.seurat, reduction = "umap")
select.cells <- CellSelector(plot = plot)

# Create a new object to fool around with
extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat <- final.pop.call.from_full.integrated.mye.seurat
Idents(extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat, cells = select.cells) <- "extra.mono"

# Check marker genes
extra.mono.markers <- FindAllMarkers(extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat, only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
extra.mono.markers <- extra.mono.markers[extra.mono.markers$cluster == "extra.mono",]
bunchOfCustomPlots(object         = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat, 
                   features       = extra.mono.markers[1:25, "gene"], 
                   assay          = "RNA", 
                   Vln.draw.names = F, 
                   Vln.width      = 25, 
                   Vln.height     = 25, 
                   Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors), 
                   ncol           = 5, 
                   name           = "final_mac_pops/extra.mono.markers")

# Display the extra cluster
customUMAP(object      = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat, 
           legend.pos  = "right", 
           file.name   = "final_mac_pops/extra.mono.UMAP.pdf", 
           pt.size     = 2,
           plot.height = 10,
           plot.width  = 17.5, 
           cols        = c("extra.mono" = "grey", M.int_refined.pop.colors))

# Get ontology
x <- extra.mono.markers[order(extra.mono.markers$gene),]
d <- duplicated(x$gene)
x <- x[!d,]
get_ontology(res = x, name = "extra.mono", outdir = "final_mac_pops/", full_GSEA = F, universe = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat, volcano.plot = F)

# Check some cell type defining genes
bunchOfCustomPlots(features       = c("CD3E", "CD4", "CD8A","FOXP3", "NCAM1", 
                                      "CD14", "CD68", "CD1C", "TREM2", "IL1B",
                                      "CD79A", "ACTA2", "CD34", "KIT", "OLR1" ), 
                   object         = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat,
                   ncol           = 5,
                   Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors), 
                   name           = "final_mac_pops/extra.mono.Cell type defining genes",
                   Vln.width      = 25, 
                   Vln.height     = 15,
                   Vln.draw.names = F
)

# Myeloid Markers
bunchOfCustomPlots(features       = c("CD200R1", "TGM2", "CD163"), 
                   object         = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat,
                   name           = "final_mac_pops/extra.mono.M2 genes",
                   Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors),  
                   Vln.draw.names = F)

bunchOfCustomPlots(features       = c("CD2","ITGAM", "PECAM1", "NCAM1", "SELL", "CSF1R"), 
                   object         = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat,
                   name           = "final_mac_pops/extra.mono.Monocyte genes", 
                   Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors), 
                   Vln.draw.names = F
)

bunchOfCustomPlots(features       = c("CD68","CD14","ITGAM","CD1C", "CLEC10A", "FCER1A", "CD2", "CD3E", "CD4"),
                   object         = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat,
                   name           = "final_mac_pops/extra.mono.Myeloid genes", 
                   ncol           = 3,
                   Vln.width      = 25, Vln.height = 15, 
                   Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors), 
                   Vln.draw.names = F
)

bunchOfCustomPlots(features       = c("CASP1","IL1B","TLR4","ABCA1","OLR1", "TREM2","CASP4","KLF4","TNF","IL18","ABCG1", "CD9"), 
                   object         = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat,
                   name           = "final_mac_pops/extra.mono.Foam vs. pro inf genes", 
                   ncol           = 3,
                   Vln.width      = 30,
                   Vln.height     = 15, Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors), 
                   Vln.draw.names = F
)

bunchOfCustomPlots(features       = c("CD14","CD68","CD86","CD1C", "FCGR2B", "CLEC10A", "FCER1A", "HLA-DQA2", "HLA-DQB2"), 
                   object         = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat,
                   name           = "final_mac_pops/extra.mono.DC genes", 
                   ncol           = 3,
                   Vln.width      = 25, 
                   Vln.height     = 15,
                   Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors),  
                   Vln.draw.names = F
)

bunchOfCustomPlots(features = c("NLRP3", "IL1B", "CASP1", "CASP4", "IL18", "PYCARD"), 
                   object   = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat,
                   name     = "final_mac_pops/extra.mono.Inflammasome genes", 
                   ncol     = 3, Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors), Vln.draw.names = F
)

bunchOfCustomPlots(features       = c("CD14", "FCGR3A", "FCGR1A", "CCR2", "PADI4", "SECISBP2L", "CD86", "ITGAX", "LYZ", "CD33", "CD36", "ITGAM"), 
                   object         = extra.mono.cluster.final.pop.call.from_full.integrated.mye.seurat,
                   name           = "final_mac_pops/extra.mono.Mono subtype genes CD14 CD16 CD64 CCR2 PADI4 SLAN CD86 CD11c LYZ CD33 CD36 CD11b", 
                   ncol           = 3,
                   Vln.width      = 25,
                   Vln.color      = c("extra.mono" = "grey", M.int_refined.pop.colors), 
                   Vln.height     = 15, 
                   Vln.draw.names = F
)

## Extra mono cluster turns out to be CD2+NKG7+PRF1+KLRD1+ NK cells, let's remove it
plot         <- DimPlot(final.pop.call.from_full.integrated.mye.seurat, reduction = "umap")

# Select all cells except the NK contaminent cluster
select.cells <- CellSelector(plot = plot)

# And subset the object
final.pop.call.from_full.integrated.mye.seurat <- subset(final.pop.call.from_full.integrated.mye.seurat, cells = select.cells)

# Update the UMAP
final.pop.call.from_full.integrated.mye.seurat <- RunUMAP(final.pop.call.from_full.integrated.mye.seurat, dims = 1:50, return.model = T)

# Draw the new UMAP
customUMAP(object     = final.pop.call.from_full.integrated.mye.seurat, 
           pt.size    = 3, 
           label      = F, 
           cols       = M.int_refined.pop.colors, 
           shuffle    = T, 
           legend.pos = "right",
           file.name  = "final_mac_pops/population UMAP plus mono clean.pdf", 
           plot.width = 15)

## Now do the total dataset
## Merge the two barely different inf pops (CD16- and S100A8+) into a fresh seurat object to keep things safe
# Check the old idents
unique(Idents(integrated.full.seurat))

# And rename
final.pop.call.integrated.full.seurat <- RenameIdents(object = integrated.full.seurat, "CD14+IL1B+SELL+CD16- Antigen-presenting Inflammatory Monocyte-derived Macrophages" = "CD14+IL1B+SELL+S100A8+ Migrating Inflammatory Monocyte-derived Macrophages")

# Check the new idents
unique(Idents(final.pop.call.integrated.full.seurat))

# Update the UMAP
final.pop.call.integrated.full.seurat <- RunUMAP(final.pop.call.integrated.full.seurat, dims = 1:50, return.model = T)

# Draw the new UMAP
customUMAP(object     = final.pop.call.integrated.full.seurat, 
           pt.size    = 3, 
           label      = F, 
           cols       = full_set.colors, 
           shuffle    = T, 
           legend.pos = "right",
           file.name  = "final_mac_pops/full set population UMAP.pdf", 
           plot.width = 20)



##=================================================================================================================
## Continue with the mac only object
## Merge the two barely different inf pops (CD16- and S100A8+) into a fresh seurat object to keep things safe
# Check the old idents
unique(Idents(integrated.mye.seurat))

# And rename
final.pop.call.integrated.mye.seurat <- RenameIdents(object = integrated.mye.seurat, "CD14+IL1B+SELL+CD16- Antigen-presenting Inflammatory Monocyte-derived Macrophages" = "CD14+IL1B+SELL+S100A8+ Migrating Inflammatory Monocyte-derived Macrophages")

# Check the new idents
unique(Idents(final.pop.call.integrated.mye.seurat))

# Update the UMAP
final.pop.call.integrated.mye.seurat <- RunUMAP(final.pop.call.integrated.mye.seurat, dims = 1:50, return.model = T)

# Draw the new UMAP
customUMAP(object     = final.pop.call.integrated.mye.seurat, 
           pt.size    = 3, 
           label      = F, 
           cols       = M.int_refined.pop.colors, 
           shuffle    = T, 
           legend.pos = "right",
           file.name  = "final_mac_pops/population UMAP.pdf", 
           plot.width = 18)

##=================================================================================================================
## Recall marker genes for the new object
# Set assay to RNA for this
DefaultAssay(final.pop.call.integrated.mye.seurat) <-"RNA"

# Rescale the data
final.pop.call.integrated.mye.seurat <- ScaleData(final.pop.call.integrated.mye.seurat, features = row.names(final.pop.call.integrated.mye.seurat))

# Call markers
final.pop.call.integrated.mye.seurat.markers <- FindAllMarkers(final.pop.call.integrated.mye.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save the top markers per cluster
sep.markers   <- final.pop.call.integrated.mye.seurat.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)
sep.markers.5 <- final.pop.call.integrated.mye.seurat.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

# plot the top9 markers per cluster
num.clus <- length(unique(Idents(final.pop.call.integrated.mye.seurat)))
for(i in levels(Idents(final.pop.call.integrated.mye.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = final.pop.call.integrated.mye.seurat,
                     name      = paste("final_mac_pops/clusters/Cluster ",i, " markers", sep = ""),
                     Vln.width = 15, Vln.height = 10, Vln.draw.names = F, feature.pt.size = 0, Vln.color = M.int_refined.pop.colors)
}

# Save markers per cluster to disk
for(i in levels(Idents(final.pop.call.integrated.mye.seurat))){
  x <- final.pop.call.integrated.mye.seurat.markers[which(final.pop.call.integrated.mye.seurat.markers$cluster == i),-6]
  x <- x[,c(6,1,2,3,4,5)]
  
  # In separate txt files
  write.table(x, file = paste("final_mac_pops/clusters/Cluster ", i, " var genes.txt", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
  
  # In one spreadsheet file
  write.xlsx(x, file = paste("final_mac_pops/clusters/Marker genes.xlsx"), sheetName = i, col.names = T, row.names = F, append = T)
}

# Make a heatmap
DoHeatmap(final.pop.call.integrated.mye.seurat, features = sep.markers.5$gene, assay = "RNA", raster = F, label = F, group.colors = M.int_refined.pop.colors)
ggsave(file = "final_mac_pops/top5.heatmap.pdf", width = 10)

## Make a heatmap of the cluster ID genes
# Extract the genes from the idents
pop.id.genes <- strsplit(levels(unique(Idents(final.pop.call.integrated.mye.seurat))), " ")
pop.id.genes <- unlist(lapply(pop.id.genes, function(x)x[1]))
pop.id.genes <- unlist(strsplit(pop.id.genes, "+", fixed = T))
pop.id.genes <- unlist(strsplit(pop.id.genes, "-", fixed = T))
pop.id.genes <- pop.id.genes[!pop.id.genes == ""]

# Fix names to symbols present in the seurat object where necessary 
pop.id.genes[which(pop.id.genes == "CD16")] <- "FCGR3A"
pop.id.genes[which(pop.id.genes == "ABCA")] <- "ABCA1"
pop.id.genes[which(pop.id.genes == "ABCG")] <- "ABCG1"

# and remove the ubiquitous CD14
pop.id.genes <- pop.id.genes[!pop.id.genes == "CD14"]

# And make the heatmap
DoHeatmap(final.pop.call.integrated.mye.seurat, features = pop.id.genes, assay = "RNA", raster = F, label = F, group.colors = M.int_refined.pop.colors)
ggsave(file = "final_mac_pops/ident_marker.heatmap.pdf", width = 10)

##==================================================================================
## Get ontologies
# All markers
for(i in levels(Idents(final.pop.call.integrated.mye.seurat))){
  x <- final.pop.call.integrated.mye.seurat.markers[which(final.pop.call.integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = i, outdir = "final_mac_pops/clusters", universe = final.pop.call.integrated.mye.seurat, full_GSEA = F, volcano.plot = T) 
}

# Top 5
for(i in levels(Idents(final.pop.call.integrated.mye.seurat))){
  x <- final.pop.call.integrated.mye.seurat.markers[which(final.pop.call.integrated.mye.seurat.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 5", i, sep = " "), outdir = "final_mac_pops/clusters", universe = final.pop.call.integrated.mye.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 5) 
}

##==================================================================================
## Plot some genes

bunchOfCustomPlots(features  = c("TREM2", "CD9", "GPNMB", "LYVE1", "MRC1", "FOLR2" ,"OLR1", "ABCA1", "ABCG1", "IL1B", "TNF", "CASP1", "CD14", "FCGR1A", "FCGR3A"), 
                   object    = final.pop.call.integrated.mye.seurat,
                   ncol      = 3,
                   name      = "final_mac_pops/LAM - resident - foam - inflammatory - general markers",
                   Vln.width = 20, Vln.height = 15, Vln.draw.names = F, Vln.color = M.int_refined.pop.colors

)

bunchOfCustomPlots(features = "PLIN2", object = final.pop.call.integrated.mye.seurat, name = "final_mac_pops/PLIN2", Vln.color = M.int_refined.pop.colors)
bunchOfCustomPlots(features = "HSPB1", object = final.pop.call.integrated.mye.seurat, name = "final_mac_pops/HSPB1", Vln.color = M.int_refined.pop.colors)


##==================================================================================
## Refine population identity curation
# Current situation
customVln(features  = pop.id.genes, 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 3,
          name      = "final_mac_pops/Curated population ID genes - violin.pdf", 
          width     = 20, height = 20, draw.names = F, Vcols = M.int_refined.pop.colors, stack = F)

# First step: Reorder
new_order.pop.id.genes <- pop.id.genes[c(1:9,
                                       11,12,10,
                                       14,15,13,
                                       16,17,18,
                                       19:27)]

customVln(features  = new_order.pop.id.genes, 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 3,
          name      = "final_mac_pops/Curated population ID genes - reordered - violin.pdf", 
          width     = 20, height = 20, draw.names = F, cols = M.int_refined.pop.colors, stack = F)

# Second step: Refine
refined.pop.id.genes <- new_order.pop.id.genes

refined.pop.id.genes[which(refined.pop.id.genes == "OLR1")]  <- "PLIN2"
refined.pop.id.genes[which(refined.pop.id.genes == "ABCG1")] <- "MMP9"
refined.pop.id.genes[which(refined.pop.id.genes == "TIMP1")] <- "PLIN2"
refined.pop.id.genes[which(refined.pop.id.genes == "NLRP3")] <- "TNF"
refined.pop.id.genes[which(refined.pop.id.genes == "ABCA1")] <- "ABCG1"

refined.pop.id.genes[c(19,22,25)] <- "ABCA1"
refined.pop.id.genes[13] <- "MRC1"
refined.pop.id.genes[15] <- "PLTP"
refined.pop.id.genes[17] <- "CD9"


customVln(features  = refined.pop.id.genes, 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 3,
          name      = "final_mac_pops/Curated population ID genes - refined - violin.pdf", 
          width     = 20, height = 20, draw.names = F, cols = M.int_refined.pop.colors, stack = F)

# Tweak it a bit
refined.pop.id.genes
refined.pop.id.genes[12] <- "PLTP"
refined.pop.id.genes[18] <- "PLTP"

customVln(features  = refined.pop.id.genes, 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 3,
          name      = "final_mac_pops/Curated population ID genes - refined - violin take 2.pdf", 
          width     = 20, height = 20, draw.names = F, cols = M.int_refined.pop.colors, stack = F)


# Lam and res markers
customVln(features  = c("CD9", "TREM2", "GPNMB", "SPP1", "FOLR2", "LYVE1", "MRC1", "SIGLEC1"), 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 4,
          name      = "final_mac_pops/LAM and resident macrophage markers.pdf", 
          width     = 23, height = 7, draw.names = F, cols = M.int_refined.pop.colors, stack = F)

# Monaco Markers
customVln(features  = c("CLEC4A", "ZEB2", "IRF4", "IRF5"), 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 2,
          name      = "final_mac_pops/Monaco markers.pdf", 
          width     = 10, height = 7, draw.names = F, cols = M.int_refined.pop.colors, stack = F)

# Fibrotic genes
customVln(features  = fib.genes, 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 4,
          name      = "final_mac_pops/fibrotic markers.pdf", 
          width     = 23, height = 7, draw.names = F, cols = M.int_refined.pop.colors, stack = F)

## subset id genes and some extra markers based on Monaco's remarks
bunchOfCustomPlots(object          = final.pop.call.integrated.mye.seurat, 
                   features        = c(unique(refined.pop.id.genes), "HMOX1", "IL10", "C1QA"),
                   assay           = "RNA", 
                   Vln.draw.names  = F, 
                   feature.pt.size = 2, 
                   Vln.width       = 20, 
                   Vln.height      = 20, Vln.color = M.int_refined.pop.colors, name = "final_mac_pops/id_genes_plus_HMOX1_IL10_C1QA")

dir.create("final_mac_pops/monaco_extra", showWarnings = F, recursive = T)
for(theGene in c(unique(refined.pop.id.genes), "HMOX1", "IL10", "C1QA")){
  customFeature(object   = final.pop.call.integrated.mye.seurat, 
                features = theGene, 
                order    = T, 
                name     = paste("final_mac_pops/monaco_extra/", theGene, " - feature plot.pdf", sep = ""), 
                pt.size  = 5)
}

FeaturePlot(object = final.pop.call.integrated.mye.seurat, features = c("PLIN2", "IL1B"), pt.size = 5, order = T, blend = T)
ggsave("final_mac_pops/monaco_extra/IL1B PLIN2 blended feature plot.pdf")

# Check CCR7
bunchOfCustomPlots(object          = final.pop.call.from_full.integrated.mye.seurat, 
                   features        = "CCR7",
                   assay           = "RNA", 
                   Vln.draw.names  = F, 
                   feature.pt.size = 2, 
                   Vln.width       = 20, 
                   Vln.height      = 20, Vln.color = M.int_refined.pop.colors, name = "final_mac_pops/CCR7")


## Run monocle analysis
## Import the seurat object
#  Convert to monocle cell data set object using SeuratWrappers
final.pop.call.from_full.integrated.mye.seurat <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = Idents(final.pop.call.from_full.integrated.mye.seurat), col.name = "final.seu.idents")
final.pop.call.integrated.mye.monocle <- as.cell_data_set(final.pop.call.from_full.integrated.mye.seurat, assay = "RNA")

## Build the single-cell trajectories
# Cluster cells with monocle
final.pop.call.integrated.mye.monocle <- cluster_cells(final.pop.call.integrated.mye.monocle, reduction_method = "UMAP", resolution = 0.01)

# Determine trajectories
final.pop.call.integrated.mye.monocle <- learn_graph(final.pop.call.integrated.mye.monocle, verbose = T, use_partition = T, learn_graph_control = list(orthogonal_proj_tip = T,
                                                                                                                     rann.k = NULL, 
                                                                                                                     prune_graph = T,
                                                                                                                     minimal_branch_len = 10,
                                                                                                                     geodesic_distance_ratio = 1/2.5,
                                                                                                                     euclidean_distance_ratio = 5))

# plot by ident
plot_cells(final.pop.call.integrated.mye.monocle, 
           label_groups_by_cluster = T, cell_stroke = F,
           label_leaves = T, 
           label_branch_points = T, 
           color_cells_by = "final.seu.idents",
           label_cell_groups = F, label_roots = T, cell_size = 3, trajectory_graph_color = "red", trajectory_graph_segment_size = 2) + 
  scale_colour_manual(breaks =  names(M.int_refined.pop.colors), values = M.int_refined.pop.colors)
ggsave(file = "final_mac_pops/monocle_trajectory.idents.pdf", width = 15, height = 10)


## Plot by pseudotime
# Set principle node as root
plot_cells(final.pop.call.integrated.mye.monocle, 
           label_principal_points = T)
ggsave(file = "final_mac_pops/monocle_trajectory_princple_points.pdf", limitsize = F)

final.pop.call.integrated.mye.monocle <- order_cells(final.pop.call.integrated.mye.monocle, root_pr_nodes = "Y_242")

plot_cells(final.pop.call.integrated.mye.monocle, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_roots = F,
           label_branch_points = FALSE, cell_size = 2, show_trajectory_graph = F)
ggsave(file = "final_mac_pops/monocle_pseudotime.pr_node_root.pdf", width = 10, height = 10)


# Set cells from classical monocyte idents to root
clas.mono.cells <- WhichCells(final.pop.call.from_full.integrated.mye.seurat, 
                              idents = unique(Idents(final.pop.call.from_full.integrated.mye.seurat))[grep("Clas", unique(Idents(final.pop.call.from_full.integrated.mye.seurat)))]
)
final.pop.call.integrated.mye.monocle <- order_cells(final.pop.call.integrated.mye.monocle, root_cells = clas.mono.cells)

plot_cells(final.pop.call.integrated.mye.monocle, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_roots = F,
           label_branch_points = FALSE, cell_size = 2, show_trajectory_graph = FALSE)
ggsave(file = "final_mac_pops/monocle_pseudotime.clas_mon_root.pdf", width = 10, height = 10)



##================================================================================================================================
##================================================================================================================================
## Run velocyto analysis
## Prep objects
# Read velocyto loom file
all.samples.loom <- ReadVelocity(file = "raw_data/merged.loom")

# Convert to Seurat object
all.samples.velocyto.seurat <- as.Seurat(x = all.samples.loom)

##================================================================================================================================
## Keep only the cells we still have in our filtered, final mye object
# Check Cell names
colnames(final.pop.call.from_full.integrated.mye.seurat)[1:5]
tail(colnames(final.pop.call.from_full.integrated.mye.seurat))
colnames(all.samples.velocyto.seurat)[1:5]
tail(all.samples.velocyto.seurat)

# Convert velocyto object cell names to the format in our 'normal' objects
velo.cell.names <- colnames(all.samples.velocyto.seurat)
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g002_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_1")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g003_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_1")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g004_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_2")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g005_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_2")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g006_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_3")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g007_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_3")

# Rename the velo cells
all.samples.velocyto.seurat <- RenameCells(all.samples.velocyto.seurat, new.names = velo.cell.names)

# And subset the object to contain only the filtered cells
all.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, cells = colnames(final.pop.call.from_full.integrated.mye.seurat))

# Same for features
DefaultAssay(final.pop.call.from_full.integrated.mye.seurat) <- "RNA"
all.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, features = row.names(final.pop.call.from_full.integrated.mye.seurat))

# And create a new mye object subsetting vice versa to contain only the cells and features present in the velo object
final.pop.call.from_full.integrated.mye.velocyto.seurat <- subset(final.pop.call.from_full.integrated.mye.seurat, cells = colnames(all.samples.velocyto.seurat))
final.pop.call.from_full.integrated.mye.velocyto.seurat <- subset(final.pop.call.from_full.integrated.mye.velocyto.seurat, features = row.names(all.samples.velocyto.seurat))


all.samples.velocyto.seurat@assays$spliced
final.pop.call.from_full.integrated.mye.velocyto.seurat@assays$RNA

##================================================================================================================================
## Transfer velocyto data to a new mye object
# Assays
final.pop.call.from_full.integrated.mye.velocyto.seurat[["spliced"]]   <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "spliced"))
final.pop.call.from_full.integrated.mye.velocyto.seurat[["unspliced"]] <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "unspliced"))
final.pop.call.from_full.integrated.mye.velocyto.seurat[["ambiguous"]] <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "ambiguous"))

# Metadata
final.pop.call.from_full.integrated.mye.velocyto.seurat <- AddMetaData(final.pop.call.from_full.integrated.mye.velocyto.seurat, metadata = all.samples.velocyto.seurat@meta.data)
head(final.pop.call.from_full.integrated.mye.velocyto.seurat@meta.data)
tail(final.pop.call.from_full.integrated.mye.velocyto.seurat@meta.data)

# Save the object
saveRDS(final.pop.call.from_full.integrated.mye.velocyto.seurat, "Seurat_Objects/final.pop.call.mye.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
#final.pop.call.from_full.integrated.mye.velocyto.seurat <- readRDS(file = "Seurat_Objects/mye.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")

##================================================================================================================================
## Clean up to save some memory
rm(all.samples.loom)
rm(all.samples.velocyto.seurat)
rm(velo.cell.names)

##================================================================================================================================
## Run velocyto natively on myeloid populations
DefaultAssay(final.pop.call.from_full.integrated.mye.velocyto.seurat) <- "spliced"

# Run velocyto scoring
final.pop.call.from_full.integrated.mye.velocyto.seurat <- RunVelocity(object = final.pop.call.from_full.integrated.mye.velocyto.seurat, deltaT = 1, kCells = 25, kGenes = 1, fit.quantile = 0.02, ncores = 3)

# Visualise
ident.colors       <- M.int_refined.pop.colors

# Get new refined idents
cell.colors        <- as.vector(Idents(final.pop.call.from_full.integrated.mye.seurat))
names(cell.colors) <- names(Idents(final.pop.call.from_full.integrated.mye.seurat))
cell.colors        <- cell.colors[names(cell.colors) %in% names(Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat))]

mono.cells         <- as.vector(Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat))[grep("SLAN", Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat))]
names(mono.cells)  <- names(Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat))[grep("SLAN", Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat))]
cell.colors        <- c(cell.colors, mono.cells)

Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat) <- cell.colors


# Set colors
for(theIdent in unique(Idents(final.pop.call.from_full.integrated.mye.velocyto.seurat))){
  cell.colors[cell.colors == theIdent] <- ident.colors[theIdent]
}

# Save the object (with updated idents)
saveRDS(final.pop.call.from_full.integrated.mye.velocyto.seurat, "Seurat_Objects/final.pop.call.mye.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")
#final.pop.call.from_full.integrated.mye.velocyto.seurat <- readRDS(file = "Seurat_Objects/mye.velo_pbmc_plaque.integrated.mapped.refined_populations.seurat.RDS")


pdf("final_mac_pops/mye.cells.velo_embeddings_n200-40.pdf",width = 15, height = 15)
show.velocity.on.embedding.cor(emb                = Embeddings(object    = final.pop.call.from_full.integrated.mye.velocyto.seurat, 
                                                               reduction = "umap"), 
                               vel                = Tool(object = final.pop.call.from_full.integrated.mye.velocyto.seurat, 
                                                         slot = "RunVelocity"),
                               n                  = 200, 
                               scale              = "sqrt", 
                               cell.colors        = ac(x     = cell.colors, alpha = 0.5), 
                               cex                = 2,
                               arrow.scale        = 1.5, 
                               show.grid.flow     = TRUE, 
                               min.grid.cell.mass = 0.5,
                               grid.n             = 40, 
                               arrow.lwd          = 2, 
                               do.par             = FALSE, 
                               cell.border.alpha  = 0,
                               n.cores            = 4,
                               nPcs               = 30)
dev.off()


##================================================================================================================================
##================================================================================================================================
## Run velocyto analysis using the mac only object
## Prep objects
# Read velocyto loom file
all.samples.loom <- ReadVelocity(file = "raw_data/merged.loom")

# Convert to Seurat object
all.samples.velocyto.seurat <- as.Seurat(x = all.samples.loom)

##================================================================================================================================
## Keep only the cells we still have in our filtered, final mye object
# Check Cell names
colnames(final.pop.call.integrated.mye.seurat)[1:5]
tail(colnames(final.pop.call.integrated.mye.seurat))
colnames(all.samples.velocyto.seurat)[1:5]
tail(all.samples.velocyto.seurat)

# Convert velocyto object cell names to the format in our 'normal' objects
velo.cell.names <- colnames(all.samples.velocyto.seurat)
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g002_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_1")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g003_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_1")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g004_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_2")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g005_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_2")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g006_ranger_out:([A-Z]+)x", replacement = "PBMC_\\1-1_3")
velo.cell.names <- gsub(x = velo.cell.names, pattern = "g007_ranger_out:([A-Z]+)x", replacement = "plaque_\\1-1_3")

# Rename the velo cells
all.samples.velocyto.seurat <- RenameCells(all.samples.velocyto.seurat, new.names = velo.cell.names)

# And subset the object to contain only the filtered cells
all.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, cells = colnames(final.pop.call.integrated.mye.seurat))

# Same for features
DefaultAssay(final.pop.call.integrated.mye.seurat) <- "RNA"
all.samples.velocyto.seurat <- subset(all.samples.velocyto.seurat, features = row.names(final.pop.call.integrated.mye.seurat))

# And create a new mye object subsetting vice versa to contain only the cells and features present in the velo object
final.pop.call.integrated.mye.velocyto.seurat <- subset(final.pop.call.integrated.mye.seurat, cells = colnames(all.samples.velocyto.seurat))
final.pop.call.integrated.mye.velocyto.seurat <- subset(final.pop.call.integrated.mye.velocyto.seurat, features = row.names(all.samples.velocyto.seurat))


all.samples.velocyto.seurat@assays$spliced
final.pop.call.integrated.mye.velocyto.seurat@assays$RNA

##================================================================================================================================
## Transfer velocyto data to a new mye object
# Assays
final.pop.call.integrated.mye.velocyto.seurat[["spliced"]]   <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "spliced"))
final.pop.call.integrated.mye.velocyto.seurat[["unspliced"]] <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "unspliced"))
final.pop.call.integrated.mye.velocyto.seurat[["ambiguous"]] <- CreateAssayObject(data = GetAssayData(all.samples.velocyto.seurat, assay = "ambiguous"))

# Metadata
final.pop.call.integrated.mye.velocyto.seurat <- AddMetaData(final.pop.call.integrated.mye.velocyto.seurat, metadata = all.samples.velocyto.seurat@meta.data)
head(final.pop.call.integrated.mye.velocyto.seurat@meta.data)
tail(final.pop.call.integrated.mye.velocyto.seurat@meta.data)

# Save the object
saveRDS(final.pop.call.integrated.mye.velocyto.seurat, "Seurat_Objects/final.pop.call.mye.velo_plaque.integrated.mapped.refined_populations.seurat.RDS")
#final.pop.call.integrated.mye.velocyto.seurat <- readRDS(file = "Seurat_Objects/mye.velo_plaque.integrated.mapped.refined_populations.seurat.RDS")

##================================================================================================================================
## Clean up to save some memory
rm(all.samples.loom)
rm(all.samples.velocyto.seurat)
rm(velo.cell.names)

##================================================================================================================================
## Run velocyto natively on myeloid populations
DefaultAssay(final.pop.call.integrated.mye.velocyto.seurat) <- "spliced"

# Run velocyto scoring
final.pop.call.integrated.mye.velocyto.seurat <- RunVelocity(object = final.pop.call.integrated.mye.velocyto.seurat, deltaT = 1, kCells = 25, kGenes = 1, fit.quantile = 0.02, ncores = 3)

# Visualise
ident.colors       <- M.int_refined.pop.colors

# Get new refined idents
cell.colors        <- as.vector(Idents(final.pop.call.integrated.mye.seurat))
names(cell.colors) <- names(Idents(final.pop.call.integrated.mye.seurat))
cell.colors        <- cell.colors[names(cell.colors) %in% names(Idents(final.pop.call.integrated.mye.velocyto.seurat))]

Idents(final.pop.call.integrated.mye.velocyto.seurat) <- cell.colors

# Set colors
for(theIdent in unique(Idents(final.pop.call.integrated.mye.velocyto.seurat))){
  cell.colors[cell.colors == theIdent] <- ident.colors[theIdent]
}

# Save the object (with updated idents)
saveRDS(final.pop.call.integrated.mye.velocyto.seurat, "Seurat_Objects/final.pop.call.mye.velo_plaque.integrated.mapped.refined_populations.seurat.RDS")
#final.pop.call.integrated.mye.velocyto.seurat <- readRDS(file = "Seurat_Objects/mye.velo_plaque.integrated.mapped.refined_populations.seurat.RDS")


pdf("final_mac_pops/mac.cells.velo_embeddings_n200-40.pdf",width = 15, height = 15)
show.velocity.on.embedding.cor(emb                = Embeddings(object    = final.pop.call.integrated.mye.velocyto.seurat, 
                                                               reduction = "umap"), 
                               vel                = Tool(object = final.pop.call.integrated.mye.velocyto.seurat, 
                                                         slot = "RunVelocity"),
                               n                  = 200, 
                               scale              = "sqrt", 
                               cell.colors        = ac(x     = cell.colors, alpha = 0.5), 
                               cex                = 2,
                               arrow.scale        = 1.5, 
                               show.grid.flow     = TRUE, 
                               min.grid.cell.mass = 0.5,
                               grid.n             = 40, 
                               arrow.lwd          = 2, 
                               do.par             = FALSE, 
                               cell.border.alpha  = 0,
                               n.cores            = 4,
                               nPcs               = 30)
dev.off()



##================================================================================================================================
##================================================================================================================================
## Check archetype markers
## Set up archetypes
archetypes <- as.vector(Idents(final.pop.call.integrated.mye.velocyto.seurat))
archetypes[grep("Foamy", archetypes)]        <- "Foamy"
archetypes[grep("Resident", archetypes)]     <- "Resident"
archetypes[grep("Lipid", archetypes)]        <- "Resident"
archetypes[grep("Monocyte-derived", archetypes)] <- "Inflammatory"

# Sanity check
levels(as.factor(archetypes))

# Add to the seurat object
final.pop.call.integrated.mye.velocyto.seurat <- AddMetaData(final.pop.call.integrated.mye.velocyto.seurat,metadata = archetypes, col.name = "archetype")

# Set up colors
archetype.colors <- c("Resident" = "#528347", "Inflammatory" = "#87482F", "Foamy" = "#9476A3", "Monocytes" = "#6495ED")

## Set up archetypes as default idents in a new object
archetype.integrated.mye.velocyto.seurat         <- AddMetaData(final.pop.call.integrated.mye.velocyto.seurat, metadata = Idents(final.pop.call.integrated.mye.velocyto.seurat), col.name = "final.pop.idents")
Idents(archetype.integrated.mye.velocyto.seurat) <- archetypes

DimPlot(archetype.integrated.mye.velocyto.seurat, cols = archetype.colors) + NoLegend()

## Calculate markers
archetype.markers <- FindAllMarkers(archetype.integrated.mye.velocyto.seurat, assay = "RNA", min.pct = 0.25, min.diff.pct = 0.25, only.pos = T)
head(archetype.markers)

# Save the top9 markers per cluster
sep.markers <- archetype.markers %>% group_by(cluster) %>% top_n(9, avg_log2FC)

# plot the top9 markers per cluster
dir.create("archetypes", showWarnings = F, recursive = T)
num.clus <- length(unique(Idents(archetype.integrated.mye.velocyto.seurat)))
for(i in levels(Idents(archetype.integrated.mye.velocyto.seurat))){
  # Make feature, violin, and dot plots
  bunchOfCustomPlots(features  = as.vector(unlist(sep.markers[sep.markers$cluster == i, "gene"])), 
                     object    = archetype.integrated.mye.velocyto.seurat,
                     name      = paste("archetypes/",i, " markers", sep = ""),
                     Vln.color = archetype.colors,
                     Vln.width = 25, Vln.height = 25)
}

# Save the UMAP for reference
customUMAP(archetype.integrated.mye.velocyto.seurat,
           pt.size = 2, 
           label = F,
           cols = archetype.colors, 
           title = "Macrophage archetypes", 
           shuffle = T, 
           legend.pos = "right", 
           seed = 666,
           file.name = "archetypes/Macrophage archetypes UMAP.pdf", 
           plot.width = 12)

# Save the markers
saveRDS(archetype.markers, file = "archetypes/Macrophage archetype markers.Rds")
saveRDS(archetype.integrated.mye.velocyto.seurat, file = "archetypes/Macrophage archetypes seurat object.Rds")
saveRDS(archetype.colors, file = "archetypes/Macrophage archetype colors.Rds")

## Make some extra plots
# Heatmap
top.markers <- archetype.markers %>% group_by(cluster) %>% top_n(15, avg_log2FC)
top.markers <- as.vector(unlist(top.markers[,"gene"]))
DoHeatmap(archetype.integrated.mye.velocyto.seurat, group.colors = archetype.colors, label = F, features = top.markers, raster = F)
ggsave("archetypes/heatmap.pdf")

# Selected markers
selected.markers <- c("PLTP", "GPNMB", "TREM2", 
                      "IFITM2", "TNFSF10", "S100A8",
                      "VEGFA", "PLIN2", "OLR1")
bunchOfCustomPlots(features  = selected.markers, 
                   object    = archetype.integrated.mye.velocyto.seurat,
                   name      = paste("archetypes/selected markers", sep = ""),
                   Vln.color = archetype.colors,
                   Vln.width = 25, Vln.height = 25)

customVln(features  = unique(refined.pop.id.genes), 
                   object    = archetype.integrated.mye.velocyto.seurat,
                   ncol = 3,
                   name      = "archetypes/refined markers.pdf",
                   col = archetype.colors,
                   width = 25, height = 25)

## Get ontologies
# Top 5
for(i in levels(Idents(archetype.integrated.mye.velocyto.seurat))){
  x <- archetype.markers[which(archetype.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 5", i, sep = " "), outdir = "archetypes", universe = archetype.integrated.mye.velocyto.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 5) 
}

# Top 10
for(i in levels(Idents(archetype.integrated.mye.velocyto.seurat))){
  x <- archetype.markers[which(archetype.markers$cluster == i),-6]
  get_ontology(x, name = paste("Top 10", i, sep = " "), outdir = "archetypes", universe = archetype.integrated.mye.velocyto.seurat, full_GSEA = F, volcano.plot = F, plot.top.n = 10) 
}

