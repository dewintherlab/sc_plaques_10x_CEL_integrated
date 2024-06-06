#======================================================================
## IPA URA
#======================================================================
## Create dir
dir.create("IPA", showWarnings = F)

## Export relevant marker lists
## Per archetype
# Write out sig markers to files
write.table(x = Mye.type.markers$Foamy,        file = "IPA/Foamy_markers.txt",        quote = F, row.names = F,col.names = F)
write.table(x = Mye.type.markers$Resident,     file = "IPA/Resident_markers.txt",     quote = F, row.names = F,col.names = F)
write.table(x = Mye.type.markers$Inflammatory, file = "IPA/Inflammatory_markers.txt", quote = F, row.names = F,col.names = F)

## Per population
final.pop.call.from_full.integrated.mye.seurat.markers <- FindAllMarkers(final.pop.call.from_full.integrated.mye.seurat, only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)

# Subset significant markers
final.pop.call.from_full.integrated.mye.seurat.sig.markers <- list()
for(i in unique(final.pop.call.from_full.integrated.mye.seurat.markers$cluster)){
  final.pop.call.from_full.integrated.mye.seurat.sig.markers[[i]] <- subset(final.pop.call.from_full.integrated.mye.seurat.markers, cluster == i & p_val_adj < 0.1)[,"gene"]
  cat(c(i, length(final.pop.call.from_full.integrated.mye.seurat.sig.markers[[i]]), "\n"))
}

# Write out sig markers to files
for(i in unique(final.pop.call.from_full.integrated.mye.seurat.markers$cluster)){
  write.table(x = final.pop.call.from_full.integrated.mye.seurat.sig.markers[[i]], file = paste("IPA/", i, ".txt", sep = ""), quote = F, row.names = F,col.names = F)
}


## Get DE genes for individual clusters -- use integrated assay so the 10X and CEL-seq smaples are better aligned
dir.create("Pseudobulk.DE/results", showWarnings = F, recursive = T)

# Remove monocyte clusters, look at the macs only
final.pop.call.from_full.integrated.mac.seurat <- subset(final.pop.call.from_full.integrated.mye.seurat, idents = as.vector(unique(Idents(final.pop.call.from_full.integrated.mye.seurat))[-grep("Monocytes",unique(Idents(final.pop.call.from_full.integrated.mye.seurat)))]))

# Prep a cluster + sample name grouping to make the pseudobulk counts over
cluster.patient <- paste(Idents(final.pop.call.from_full.integrated.mac.seurat), final.pop.call.from_full.integrated.mac.seurat$Patient, sep = ".")
cluster.patient <- gsub(pattern = " ", replacement = "_", x = cluster.patient)
cluster.patient <- gsub(pattern = "+", replacement = "hi_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "-", replacement = "lo_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Nonlo", replacement = "Non", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Monocytelo", replacement = "Monocyte_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Residentlo", replacement = "Resident_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "__", replacement = "_", x = cluster.patient, fixed = T)

final.pop.call.from_full.integrated.mac.seurat <- AddMetaData(final.pop.call.from_full.integrated.mac.seurat, metadata = cluster.patient, col.name = "clus.sample")

# Make pseudobulk
final.pop.call.from_full.integrated.mac.avgexp <- AverageExpression(object = final.pop.call.from_full.integrated.mac.seurat, group.by = "clus.sample", assays = "integrated")$integrated

# Make integers (to simulate readcounts)
final.pop.call.from_full.integrated.mac.avgexp <- ceiling(final.pop.call.from_full.integrated.mac.avgexp)

# Make metadata table
DE.metadata <- data.frame("Cluster" = t(as.data.frame(strsplit(colnames(final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,1], 
                          "Patient" = t(as.data.frame(strsplit(colnames(final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,2], 
                          row.names = colnames(final.pop.call.from_full.integrated.mac.avgexp))
colnames(DE.metadata) <- c("Cluster", "Patient")
DE.metadata$Cluster <- as.factor(DE.metadata$Cluster)
DE.metadata$Patient <- as.factor(DE.metadata$Patient)
head(DE.metadata)

# setup colors
pseudo.DE.colors <- M.int_refined.pop.colors
names(pseudo.DE.colors) <- gsub(pattern = " ", replacement = "_", x = names(pseudo.DE.colors))
names(pseudo.DE.colors) <- gsub(pattern = "+", replacement = "hi_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "-", replacement = "lo_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Nonlo", replacement = "Non", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Monocytelo", replacement = "Monocyte_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Residentlo", replacement = "Resident_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "__", replacement = "_", x = names(pseudo.DE.colors), fixed = T)
pseudo.DE.colors <- pseudo.DE.colors[-c(1,2,3)]


## Run DESeq2 for DE analyses
# Get the data ready
dds <- DESeqDataSetFromMatrix(countData = final.pop.call.from_full.integrated.mac.avgexp, colData = DE.metadata, design = ~ Patient + Cluster)

# Remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds, parallel = T)
vst <- varianceStabilizingTransformation(dds)
dds
vst

## Visualise general clustering
topVarGenes <- head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = DE.metadata,
         annotation_colors = list("Cluster" = pseudo.DE.colors),
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "Pseudobulk.DE/top1000_variable_genes.kmeans_heatmap.pdf", width = 20, height = 20
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$Patient <- DE.metadata$Patient 
mat.pca$Cluster <- DE.metadata$Cluster
percentVar <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Cluster)) +
  geom_point(size=3) +
  scale_color_manual(values = pseudo.DE.colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE/top1000_variable_genes.PCA_plot_Cluster.pdf", width = 20, height = 10)

ggplot(mat.pca, aes(PC1, PC2, color = Patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE/top1000_variable_genes.PCA_plot_Patient.pdf")

resultsNames(dds)

res.int <- list()
for(i in levels(DE.metadata$Cluster)){
  for(j in unique(DE.metadata$Cluster)){
    if(i != j){ # Don't compare to self
      if(! paste(j, "_over_", i, sep = "") %in% names(res.int)){ # We don't need to do the reciprocals
        res.int[[paste(i, "_over_", j, sep = "")]] <- results(dds, contrast = c("Cluster", i, j), parallel = T)
        write.table(row.names(subset(res.int[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0)), file = paste("Pseudobulk.DE/results/", i, "_over_", j, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        write.table(row.names(subset(res.int[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0)), file = paste("Pseudobulk.DE/results/", j, "_over_", i, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        get_ontology_DE(res = subset(res.int[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0), name = paste(i, "_over_", j, sep = ""), outdir = "Pseudobulk.DE/results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
        get_ontology_DE(res = subset(res.int[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0), name = paste(j, "_over_", i, sep = ""), outdir = "Pseudobulk.DE/results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
      }
    }
  }
}

# Order the res tables on significance
for(i in names(res.int)){
  res.int[[i]] <- res.int[[i]][order(res.int[[i]]$padj),]
}

# Plot the top 50 DE genes
for(i in names(res.int)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.int[[i]])[1:50],])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE/results/",i, ".top50_DE.scaled_heatmap.pdf", sep = "")
           )
}

# Plot the sig DE genes
for(i in names(res.int)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.int[[i]], padj < 0.1)
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors), 
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE/results/",i, ".sig_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the top 50 DE genes - in the affected clusters only
for(i in names(res.int)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.int[[i]])[1:50],])
  a      <- strsplit(i, "_over_")[[1]][[1]]
  b      <- strsplit(i, "_over_")[[1]][[2]]
  keep   <- grep(a, colnames(de.mat))
  keep   <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE/results/",i, ".top50_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.int)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.int[[i]], padj < 0.1)
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  a              <- strsplit(i, "_over_")[[1]][[1]]
  b              <- strsplit(i, "_over_")[[1]][[2]]
  keep           <- grep(a, colnames(de.mat))
  keep           <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE/results/",i, ".sig_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}


##===========================================================================
## Get DE genes for archetypes
dir.create("Pseudobulk.DE/archetype_results", showWarnings = F, recursive = T)

## Set up archetypes
archetypes <- as.vector(Idents(final.pop.call.from_full.integrated.mye.seurat))
archetypes[grep("Foamy", archetypes)]        <- "Foamy"
archetypes[grep("Resident", archetypes)]     <- "Resident"
archetypes[grep("Lipid", archetypes)]        <- "Resident"
archetypes[grep("Monocyte-derived", archetypes)] <- "Inflammatory"

# Sanity check
levels(as.factor(archetypes))

# Add to the seurat object
final.pop.call.from_full.integrated.mye.seurat <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = archetypes, col.name = "archetype")

# Set up colors
archetype.colors <- c("Resident" = "#528347", "Inflammatory" = "#87482F", "Foamy" = "#9476A3", "Monocytes" = "#6495ED")

## Set up archetypes as default idents in a new object
archetype.integrated.mye.velocyto.seurat         <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = Idents(final.pop.call.from_full.integrated.mye.seurat), col.name = "final.pop.idents")
Idents(archetype.integrated.mye.velocyto.seurat) <- archetypes

# Remove monocyte clusters, look at the macs only
archetype.final.pop.call.from_full.integrated.mac.seurat <- subset(archetype.integrated.mye.velocyto.seurat, idents = as.vector(unique(Idents(archetype.integrated.mye.velocyto.seurat))[-grep("Monocytes",unique(Idents(archetype.integrated.mye.velocyto.seurat)))]))
unique(Idents(archetype.final.pop.call.from_full.integrated.mac.seurat))

# Prep a cluster + sample name grouping to make the pseudobulk counts over
cluster.patient <- paste(Idents(archetype.final.pop.call.from_full.integrated.mac.seurat), archetype.final.pop.call.from_full.integrated.mac.seurat$Patient, sep = ".")

archetype.final.pop.call.from_full.integrated.mac.seurat <- AddMetaData(archetype.final.pop.call.from_full.integrated.mac.seurat, metadata = cluster.patient, col.name = "clus.sample")

# Make pseudobulk
archetype.final.pop.call.from_full.integrated.mac.avgexp <- AverageExpression(object = archetype.final.pop.call.from_full.integrated.mac.seurat, group.by = "clus.sample", assays = "integrated")$integrated

# Make integers (to simulate readcounts)
archetype.final.pop.call.from_full.integrated.mac.avgexp <- ceiling(archetype.final.pop.call.from_full.integrated.mac.avgexp)

# Make metadata table
DE.metadata <- data.frame("Cluster" = t(as.data.frame(strsplit(colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,1], 
                          "Patient" = t(as.data.frame(strsplit(colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,2], 
                          row.names = colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp))
colnames(DE.metadata) <- c("Cluster", "Patient")
DE.metadata$Cluster <- as.factor(DE.metadata$Cluster)
DE.metadata$Patient <- as.factor(DE.metadata$Patient)
head(DE.metadata)

# setup colors
pseudo.DE.colors <- archetype.colors

## Run DESeq2 for DE analyses
# Get the data ready
dds <- DESeqDataSetFromMatrix(countData = archetype.final.pop.call.from_full.integrated.mac.avgexp, colData = DE.metadata, design = ~ Patient + Cluster)

# Remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds, parallel = T)
vst <- varianceStabilizingTransformation(dds)
dds
vst

## Visualise general clustering
topVarGenes <- head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = DE.metadata,
         annotation_colors = list("Cluster" = pseudo.DE.colors),
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "Pseudobulk.DE/archetype_results/top1000_variable_genes.kmeans_heatmap.pdf", width = 20, height = 20
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$Patient <- DE.metadata$Patient 
mat.pca$Cluster <- DE.metadata$Cluster
percentVar <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Cluster)) +
  geom_point(size=3) +
  scale_color_manual(values = pseudo.DE.colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE/archetype_results/top1000_variable_genes.PCA_plot_Cluster.pdf", width = 20, height = 10)

ggplot(mat.pca, aes(PC1, PC2, color = Patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE/archetype_results/top1000_variable_genes.PCA_plot_Patient.pdf")

resultsNames(dds)

res.arch.int <- list()
for(i in levels(DE.metadata$Cluster)){
  for(j in unique(DE.metadata$Cluster)){
    if(i != j){ # Don't compare to self
      if(! paste(j, "_over_", i, sep = "") %in% names(res.arch.int)){ # We don't need to do the reciprocals
        cat(paste(i, "_over_", j, "\n",sep = ""))
        res.arch.int[[paste(i, "_over_", j, sep = "")]] <- results(dds, contrast = c("Cluster", i, j), parallel = T)
        write.table(row.names(subset(res.arch.int[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0)), file = paste("Pseudobulk.DE/archetype_results/", i, "_over_", j, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        write.table(row.names(subset(res.arch.int[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0)), file = paste("Pseudobulk.DE/archetype_results/", j, "_over_", i, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        get_ontology_DE(res = subset(res.arch.int[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0), name = paste(i, "_over_", j, sep = ""), outdir = "Pseudobulk.DE/archetype_results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
        get_ontology_DE(res = subset(res.arch.int[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0), name = paste(j, "_over_", i, sep = ""), outdir = "Pseudobulk.DE/archetype_results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
      }
    }
  }
}

# Order the res tables on significance
for(i in names(res.arch.int)){
  res.arch.int[[i]] <- res.arch.int[[i]][order(res.arch.int[[i]]$padj),]
}

# Plot the top 50 DE genes
for(i in names(res.arch.int)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.arch.int[[i]])[1:50],])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE/archetype_results/",i, ".top50_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.arch.int)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.arch.int[[i]], padj < 0.1)
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors), 
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE/archetype_results/",i, ".sig_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the top 50 DE genes - in the affected clusters only
for(i in names(res.arch.int)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.arch.int[[i]])[1:50],])
  a      <- strsplit(i, "_over_")[[1]][[1]]
  b      <- strsplit(i, "_over_")[[1]][[2]]
  keep   <- grep(a, colnames(de.mat))
  keep   <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE/archetype_results/",i, ".top50_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.arch.int)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.arch.int[[i]], padj < 0.1)
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  a              <- strsplit(i, "_over_")[[1]][[1]]
  b              <- strsplit(i, "_over_")[[1]][[2]]
  keep           <- grep(a, colnames(de.mat))
  keep           <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE/archetype_results/",i, ".sig_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}


## Get DE genes for individual clusters -- Use RNA assay
dir.create("Pseudobulk.DE.RNA/results", showWarnings = F, recursive = T)

# Remove monocyte clusters, look at the macs only
final.pop.call.from_full.integrated.mac.seurat <- subset(final.pop.call.from_full.integrated.mye.seurat, idents = as.vector(unique(Idents(final.pop.call.from_full.integrated.mye.seurat))[-grep("Monocytes",unique(Idents(final.pop.call.from_full.integrated.mye.seurat)))]))

# Prep a cluster + sample name grouping to make the pseudobulk counts over
cluster.patient <- paste(Idents(final.pop.call.from_full.integrated.mac.seurat), final.pop.call.from_full.integrated.mac.seurat$Patient, sep = ".")
cluster.patient <- gsub(pattern = " ", replacement = "_", x = cluster.patient)
cluster.patient <- gsub(pattern = "+", replacement = "hi_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "-", replacement = "lo_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Nonlo", replacement = "Non", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Monocytelo", replacement = "Monocyte_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Residentlo", replacement = "Resident_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "__", replacement = "_", x = cluster.patient, fixed = T)

final.pop.call.from_full.integrated.mac.seurat <- AddMetaData(final.pop.call.from_full.integrated.mac.seurat, metadata = cluster.patient, col.name = "clus.sample")

# Make pseudobulk
final.pop.call.from_full.integrated.mac.avgexp <- AverageExpression(object = final.pop.call.from_full.integrated.mac.seurat, group.by = "clus.sample", assays = "RNA")$RNA

# Make integers (to simulate readcounts)
final.pop.call.from_full.integrated.mac.avgexp <- ceiling(final.pop.call.from_full.integrated.mac.avgexp)

# Make metadata table
DE.metadata <- data.frame("Cluster" = t(as.data.frame(strsplit(colnames(final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,1], 
                          "Patient" = t(as.data.frame(strsplit(colnames(final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,2], 
                          row.names = colnames(final.pop.call.from_full.integrated.mac.avgexp))
colnames(DE.metadata) <- c("Cluster", "Patient")
DE.metadata$Cluster <- as.factor(DE.metadata$Cluster)
DE.metadata$Patient <- as.factor(DE.metadata$Patient)
head(DE.metadata)

# setup colors
pseudo.DE.colors <- M.int_refined.pop.colors
names(pseudo.DE.colors) <- gsub(pattern = " ", replacement = "_", x = names(pseudo.DE.colors))
names(pseudo.DE.colors) <- gsub(pattern = "+", replacement = "hi_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "-", replacement = "lo_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Nonlo", replacement = "Non", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Monocytelo", replacement = "Monocyte_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Residentlo", replacement = "Resident_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "__", replacement = "_", x = names(pseudo.DE.colors), fixed = T)
pseudo.DE.colors <- pseudo.DE.colors[-c(1,2,3)]


## Run DESeq2 for DE analyses
# Get the data ready
dds <- DESeqDataSetFromMatrix(countData = final.pop.call.from_full.integrated.mac.avgexp, colData = DE.metadata, design = ~ Patient + Cluster)

# Remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds, parallel = T)
vst <- varianceStabilizingTransformation(dds)
dds
vst

## Visualise general clustering
topVarGenes <- head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = DE.metadata,
         annotation_colors = list("Cluster" = pseudo.DE.colors),
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "Pseudobulk.DE.RNA/top1000_variable_genes.kmeans_heatmap.pdf", width = 20, height = 20
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$Patient <- DE.metadata$Patient 
mat.pca$Cluster <- DE.metadata$Cluster
percentVar <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Cluster)) +
  geom_point(size=3) +
  scale_color_manual(values = pseudo.DE.colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.RNA/top1000_variable_genes.PCA_plot_Cluster.pdf", width = 20, height = 10)

ggplot(mat.pca, aes(PC1, PC2, color = Patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.RNA/top1000_variable_genes.PCA_plot_Patient.pdf")

resultsNames(dds)

res.RNA <- list()
for(i in levels(DE.metadata$Cluster)){
  for(j in unique(DE.metadata$Cluster)){
    if(i != j){ # Don't compare to self
      if(! paste(j, "_over_", i, sep = "") %in% names(res.RNA)){ # We don't need to do the reciprocals
        cat(paste(i, "_over_", j, "\n",sep = ""))
        res.RNA[[paste(i, "_over_", j, sep = "")]] <- results(dds, contrast = c("Cluster", i, j), parallel = T)
        write.table(row.names(subset(res.RNA[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0)), file = paste("Pseudobulk.DE.RNA/results/", i, "_over_", j, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        write.table(row.names(subset(res.RNA[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0)), file = paste("Pseudobulk.DE.RNA/results/", j, "_over_", i, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        get_ontology_DE(res = subset(res.RNA[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0), name = paste(i, "_over_", j, sep = ""), outdir = "Pseudobulk.DE.RNA/results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
        get_ontology_DE(res = subset(res.RNA[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0), name = paste(j, "_over_", i, sep = ""), outdir = "Pseudobulk.DE.RNA/results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
      }
    }
  }
}

# Order the res tables on significance
for(i in names(res.RNA)){
  res.RNA[[i]] <- res.RNA[[i]][order(res.RNA[[i]]$padj),]
}

# Plot the top 50 DE genes
for(i in names(res.RNA)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.RNA[[i]])[1:50],])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.RNA/results/",i, ".top50_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.RNA)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.RNA[[i]], padj < 0.1)
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors), 
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.RNA/results/",i, ".sig_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the top 50 DE genes - in the affected clusters only
for(i in names(res.RNA)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.RNA[[i]])[1:50],])
  a      <- strsplit(i, "_over_")[[1]][[1]]
  b      <- strsplit(i, "_over_")[[1]][[2]]
  keep   <- grep(a, colnames(de.mat))
  keep   <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.RNA/results/",i, ".top50_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.RNA)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.RNA[[i]], padj < 0.1)
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  a              <- strsplit(i, "_over_")[[1]][[1]]
  b              <- strsplit(i, "_over_")[[1]][[2]]
  keep           <- grep(a, colnames(de.mat))
  keep           <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.RNA/results/",i, ".sig_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}


##===========================================================================
## Get DE genes for archetypes
dir.create("Pseudobulk.DE.RNA/archetype_results", showWarnings = F, recursive = T)

## Set up archetypes
archetypes <- as.vector(Idents(final.pop.call.from_full.integrated.mye.seurat))
archetypes[grep("Foamy", archetypes)]        <- "Foamy"
archetypes[grep("Resident", archetypes)]     <- "Resident"
archetypes[grep("Lipid", archetypes)]        <- "Resident"
archetypes[grep("Monocyte-derived", archetypes)] <- "Inflammatory"

# Sanity check
levels(as.factor(archetypes))

# Add to the seurat object
final.pop.call.from_full.integrated.mye.seurat <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = archetypes, col.name = "archetype")

# Set up colors
archetype.colors <- c("Resident" = "#528347", "Inflammatory" = "#87482F", "Foamy" = "#9476A3", "Monocytes" = "#6495ED")

## Set up archetypes as default idents in a new object
archetype.integrated.mye.velocyto.seurat         <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = Idents(final.pop.call.from_full.integrated.mye.seurat), col.name = "final.pop.idents")
Idents(archetype.integrated.mye.velocyto.seurat) <- archetypes

# Remove monocyte clusters, look at the macs only
archetype.final.pop.call.from_full.integrated.mac.seurat <- subset(archetype.integrated.mye.velocyto.seurat, idents = as.vector(unique(Idents(archetype.integrated.mye.velocyto.seurat))[-grep("Monocytes",unique(Idents(archetype.integrated.mye.velocyto.seurat)))]))
unique(Idents(archetype.final.pop.call.from_full.integrated.mac.seurat))

# Prep a cluster + sample name grouping to make the pseudobulk counts over
cluster.patient <- paste(Idents(archetype.final.pop.call.from_full.integrated.mac.seurat), archetype.final.pop.call.from_full.integrated.mac.seurat$Patient, sep = ".")

archetype.final.pop.call.from_full.integrated.mac.seurat <- AddMetaData(archetype.final.pop.call.from_full.integrated.mac.seurat, metadata = cluster.patient, col.name = "clus.sample")

# Make pseudobulk
archetype.final.pop.call.from_full.integrated.mac.avgexp <- AverageExpression(object = archetype.final.pop.call.from_full.integrated.mac.seurat, group.by = "clus.sample", assays = "RNA")$RNA

# Make integers (to simulate readcounts)
archetype.final.pop.call.from_full.integrated.mac.avgexp <- ceiling(archetype.final.pop.call.from_full.integrated.mac.avgexp)

# Make metadata table
DE.metadata <- data.frame("Cluster" = t(as.data.frame(strsplit(colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,1], 
                          "Patient" = t(as.data.frame(strsplit(colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,2], 
                          row.names = colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp))
colnames(DE.metadata) <- c("Cluster", "Patient")
DE.metadata$Cluster <- as.factor(DE.metadata$Cluster)
DE.metadata$Patient <- as.factor(DE.metadata$Patient)
head(DE.metadata)

# setup colors
pseudo.DE.colors <- archetype.colors

## Run DESeq2 for DE analyses
# Get the data ready
dds <- DESeqDataSetFromMatrix(countData = archetype.final.pop.call.from_full.integrated.mac.avgexp, colData = DE.metadata, design = ~ Patient + Cluster)

# Remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds, parallel = T)
vst <- varianceStabilizingTransformation(dds)
dds
vst

## Visualise general clustering
topVarGenes <- head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = DE.metadata,
         annotation_colors = list("Cluster" = pseudo.DE.colors),
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "Pseudobulk.DE.RNA/archetype_results/top1000_variable_genes.kmeans_heatmap.pdf", width = 20, height = 20
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$Patient <- DE.metadata$Patient 
mat.pca$Cluster <- DE.metadata$Cluster
percentVar <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Cluster)) +
  geom_point(size=3) +
  scale_color_manual(values = pseudo.DE.colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.RNA/archetype_results/top1000_variable_genes.PCA_plot_Cluster.pdf", width = 20, height = 10)

ggplot(mat.pca, aes(PC1, PC2, color = Patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.RNA/archetype_results/top1000_variable_genes.PCA_plot_Patient.pdf")

resultsNames(dds)

res.RNA.arch <- list()
for(i in levels(DE.metadata$Cluster)){
  for(j in unique(DE.metadata$Cluster)){
    if(i != j){ # Don't compare to self
      if(! paste(j, "_over_", i, sep = "") %in% names(res.RNA.arch)){ # We don't need to do the reciprocals
        cat(paste(i, "_over_", j, "\n",sep = ""))
        res.RNA.arch[[paste(i, "_over_", j, sep = "")]] <- results(dds, contrast = c("Cluster", i, j), parallel = T)
        write.table(row.names(subset(res.RNA.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0)), file = paste("Pseudobulk.DE.RNA/archetype_results/", i, "_over_", j, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        write.table(row.names(subset(res.RNA.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0)), file = paste("Pseudobulk.DE.RNA/archetype_results/", j, "_over_", i, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        get_ontology_DE(res = subset(res.RNA.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0), name = paste(i, "_over_", j, sep = ""), outdir = "Pseudobulk.DE.RNA/archetype_results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
        get_ontology_DE(res = subset(res.RNA.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0), name = paste(j, "_over_", i, sep = ""), outdir = "Pseudobulk.DE.RNA/archetype_results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
      }
    }
  }
}

# Order the res tables on significance
for(i in names(res.RNA.arch)){
  res.RNA.arch[[i]] <- res.RNA.arch[[i]][order(res.RNA.arch[[i]]$padj),]
}

# Plot the top 50 DE genes
for(i in names(res.RNA.arch)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.RNA.arch[[i]])[1:50],])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.RNA/archetype_results/",i, ".top50_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.RNA.arch)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.RNA.arch[[i]], padj < 0.1)
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors), 
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.RNA/archetype_results/",i, ".sig_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the top 50 DE genes - in the affected clusters only
for(i in names(res.RNA.arch)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.RNA.arch[[i]])[1:50],])
  a      <- strsplit(i, "_over_")[[1]][[1]]
  b      <- strsplit(i, "_over_")[[1]][[2]]
  keep   <- grep(a, colnames(de.mat))
  keep   <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.RNA/archetype_results/",i, ".top50_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.RNA.arch)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.RNA.arch[[i]], padj < 0.1)
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  a              <- strsplit(i, "_over_")[[1]][[1]]
  b              <- strsplit(i, "_over_")[[1]][[2]]
  keep           <- grep(a, colnames(de.mat))
  keep           <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.RNA/archetype_results/",i, ".sig_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}


## Get DE genes for individual clusters -- Use CEL-seq only RNA assay
dir.create("Pseudobulk.DE.CEL.RNA/results", showWarnings = F, recursive = T)

# Remove monocyte clusters, look at the macs only
final.pop.call.from_full.integrated.mac.seurat <- subset(final.pop.call.from_full.integrated.mye.seurat, idents = as.vector(unique(Idents(final.pop.call.from_full.integrated.mye.seurat))[-grep("Monocytes",unique(Idents(final.pop.call.from_full.integrated.mye.seurat)))]))

# Remove 10X patients
final.pop.call.from_full.integrated.mac.CEL_only.seurat <- subset(final.pop.call.from_full.integrated.mac.seurat, subset = Patient != "P1" & Patient != "P2" & Patient != "P3")
unique(final.pop.call.from_full.integrated.mac.CEL_only.seurat$Patient)

# Prep a cluster + sample name grouping to make the pseudobulk counts over
cluster.patient <- paste(Idents(final.pop.call.from_full.integrated.mac.CEL_only.seurat), final.pop.call.from_full.integrated.mac.CEL_only.seurat$Patient, sep = ".")
cluster.patient <- gsub(pattern = " ", replacement = "_", x = cluster.patient)
cluster.patient <- gsub(pattern = "+", replacement = "hi_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "-", replacement = "lo_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Nonlo", replacement = "Non", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Monocytelo", replacement = "Monocyte_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Residentlo", replacement = "Resident_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "__", replacement = "_", x = cluster.patient, fixed = T)

final.pop.call.from_full.integrated.mac.CEL_only.seurat <- AddMetaData(final.pop.call.from_full.integrated.mac.CEL_only.seurat, metadata = cluster.patient, col.name = "clus.sample")

# Make pseudobulk
final.pop.call.from_full.integrated.mac.avgexp <- AverageExpression(object = final.pop.call.from_full.integrated.mac.CEL_only.seurat, group.by = "clus.sample", assays = "RNA")$RNA

# Make integers (to simulate readcounts)
final.pop.call.from_full.integrated.mac.avgexp <- ceiling(final.pop.call.from_full.integrated.mac.avgexp)

# Make metadata table
DE.metadata <- data.frame("Cluster" = t(as.data.frame(strsplit(colnames(final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,1], 
                          "Patient" = t(as.data.frame(strsplit(colnames(final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,2], 
                          row.names = colnames(final.pop.call.from_full.integrated.mac.avgexp))
colnames(DE.metadata) <- c("Cluster", "Patient")
DE.metadata$Cluster <- as.factor(DE.metadata$Cluster)
DE.metadata$Patient <- as.factor(DE.metadata$Patient)
head(DE.metadata)

# setup colors
pseudo.DE.colors <- M.int_refined.pop.colors
names(pseudo.DE.colors) <- gsub(pattern = " ", replacement = "_", x = names(pseudo.DE.colors))
names(pseudo.DE.colors) <- gsub(pattern = "+", replacement = "hi_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "-", replacement = "lo_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Nonlo", replacement = "Non", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Monocytelo", replacement = "Monocyte_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Residentlo", replacement = "Resident_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "__", replacement = "_", x = names(pseudo.DE.colors), fixed = T)
pseudo.DE.colors <- pseudo.DE.colors[-c(1,2,3)]


## Run DESeq2 for DE analyses
# Get the data ready
dds <- DESeqDataSetFromMatrix(countData = final.pop.call.from_full.integrated.mac.avgexp, colData = DE.metadata, design = ~ Patient + Cluster)

# Remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds, parallel = T)
vst <- varianceStabilizingTransformation(dds)
dds
vst

## Visualise general clustering
topVarGenes <- head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = DE.metadata,
         annotation_colors = list("Cluster" = pseudo.DE.colors),
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "Pseudobulk.DE.CEL.RNA/top1000_variable_genes.kmeans_heatmap.pdf", width = 20, height = 20
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$Patient <- DE.metadata$Patient 
mat.pca$Cluster <- DE.metadata$Cluster
percentVar <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Cluster)) +
  geom_point(size=3) +
  scale_color_manual(values = pseudo.DE.colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.CEL.RNA/top1000_variable_genes.PCA_plot_Cluster.pdf", width = 20, height = 10)

ggplot(mat.pca, aes(PC1, PC2, color = Patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.CEL.RNA/top1000_variable_genes.PCA_plot_Patient.pdf")

resultsNames(dds)

res.CEL <- list()
for(i in levels(DE.metadata$Cluster)){
  for(j in unique(DE.metadata$Cluster)){
    if(i != j){ # Don't compare to self
      if(! paste(j, "_over_", i, sep = "") %in% names(res.CEL)){ # We don't need to do the reciprocals
        cat(paste(i, "_over_", j, "\n",sep = ""))
        res.CEL[[paste(i, "_over_", j, sep = "")]] <- results(dds, contrast = c("Cluster", i, j), parallel = T)
        write.table(row.names(subset(res.CEL[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0)), file = paste("Pseudobulk.DE.CEL.RNA/results/", i, "_over_", j, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        write.table(row.names(subset(res.CEL[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0)), file = paste("Pseudobulk.DE.CEL.RNA/results/", j, "_over_", i, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        get_ontology_DE(res = subset(res.CEL[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0), name = paste(i, "_over_", j, sep = ""), outdir = "Pseudobulk.DE.CEL.RNA/results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
        get_ontology_DE(res = subset(res.CEL[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0), name = paste(j, "_over_", i, sep = ""), outdir = "Pseudobulk.DE.CEL.RNA/results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
      }
    }
  }
}

# Order the res tables on significance
for(i in names(res.CEL)){
  res.CEL[[i]] <- res.CEL[[i]][order(res.CEL[[i]]$padj),]
}

# Plot the top 50 DE genes
for(i in names(res.CEL)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.CEL[[i]])[1:50],])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.CEL.RNA/results/",i, ".top50_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.CEL)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.CEL[[i]], padj < 0.1)
  if(nrow(fc_ordered_res) < 1){next}
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors), 
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.CEL.RNA/results/",i, ".sig_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the top 50 DE genes - in the affected clusters only
for(i in names(res.CEL)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.CEL[[i]])[1:50],])
  a      <- strsplit(i, "_over_")[[1]][[1]]
  b      <- strsplit(i, "_over_")[[1]][[2]]
  keep   <- grep(a, colnames(de.mat))
  keep   <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.CEL.RNA/results/",i, ".top50_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.CEL)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.CEL[[i]], padj < 0.1)
  if(nrow(fc_ordered_res) < 1){next}
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  a              <- strsplit(i, "_over_")[[1]][[1]]
  b              <- strsplit(i, "_over_")[[1]][[2]]
  keep           <- grep(a, colnames(de.mat))
  keep           <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.CEL.RNA/results/",i, ".sig_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}


##===========================================================================
## Get DE genes for archetypes
dir.create("Pseudobulk.DE.CEL.RNA/archetype_results", showWarnings = F, recursive = T)

## Set up archetypes
archetypes <- as.vector(Idents(final.pop.call.from_full.integrated.mye.seurat))
archetypes[grep("Foamy", archetypes)]        <- "Foamy"
archetypes[grep("Resident", archetypes)]     <- "Resident"
archetypes[grep("Lipid", archetypes)]        <- "Resident"
archetypes[grep("Monocyte-derived", archetypes)] <- "Inflammatory"

# Sanity check
levels(as.factor(archetypes))

# Add to the seurat object
final.pop.call.from_full.integrated.mye.seurat <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = archetypes, col.name = "archetype")

# Set up colors
archetype.colors <- c("Resident" = "#528347", "Inflammatory" = "#87482F", "Foamy" = "#9476A3", "Monocytes" = "#6495ED")

## Set up archetypes as default idents in a new object
archetype.integrated.mye.velocyto.seurat         <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = Idents(final.pop.call.from_full.integrated.mye.seurat), col.name = "final.pop.idents")
Idents(archetype.integrated.mye.velocyto.seurat) <- archetypes

# Remove monocyte clusters, look at the macs only
archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat <- subset(archetype.integrated.mye.velocyto.seurat, idents = as.vector(unique(Idents(archetype.integrated.mye.velocyto.seurat))[-grep("Monocytes",unique(Idents(archetype.integrated.mye.velocyto.seurat)))]))
unique(Idents(archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat))

# Remove 10X patients
archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat <- subset(archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat, subset = Patient != "P1" & Patient != "P2" & Patient != "P3")
unique(archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat$Patient)

# Prep a cluster + sample name grouping to make the pseudobulk counts over
cluster.patient <- paste(Idents(archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat), archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat$Patient, sep = ".")

archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat <- AddMetaData(archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat, metadata = cluster.patient, col.name = "clus.sample")

# Make pseudobulk
archetype.final.pop.call.from_full.integrated.mac.avgexp <- AverageExpression(object = archetype.final.pop.call.from_full.integrated.mac.CEL_only.seurat, group.by = "clus.sample", assays = "RNA")$RNA

# Make integers (to simulate readcounts)
archetype.final.pop.call.from_full.integrated.mac.avgexp <- ceiling(archetype.final.pop.call.from_full.integrated.mac.avgexp)

# Make metadata table
DE.metadata <- data.frame("Cluster" = t(as.data.frame(strsplit(colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,1], 
                          "Patient" = t(as.data.frame(strsplit(colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,2], 
                          row.names = colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp))
colnames(DE.metadata) <- c("Cluster", "Patient")
DE.metadata$Cluster <- as.factor(DE.metadata$Cluster)
DE.metadata$Patient <- as.factor(DE.metadata$Patient)
head(DE.metadata)

# setup colors
pseudo.DE.colors <- archetype.colors

## Run DESeq2 for DE analyses
# Get the data ready
dds <- DESeqDataSetFromMatrix(countData = archetype.final.pop.call.from_full.integrated.mac.avgexp, colData = DE.metadata, design = ~ Patient + Cluster)

# Remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds, parallel = T)
vst <- varianceStabilizingTransformation(dds)
dds
vst

## Visualise general clustering
topVarGenes <- head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = DE.metadata,
         annotation_colors = list("Cluster" = pseudo.DE.colors),
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "Pseudobulk.DE.CEL.RNA/archetype_results/top1000_variable_genes.kmeans_heatmap.pdf", width = 20, height = 20
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$Patient <- DE.metadata$Patient 
mat.pca$Cluster <- DE.metadata$Cluster
percentVar <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Cluster)) +
  geom_point(size=3) +
  scale_color_manual(values = pseudo.DE.colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.CEL.RNA/archetype_results/top1000_variable_genes.PCA_plot_Cluster.pdf", width = 20, height = 10)

ggplot(mat.pca, aes(PC1, PC2, color = Patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.CEL.RNA/archetype_results/top1000_variable_genes.PCA_plot_Patient.pdf")

resultsNames(dds)

res.CEL.arch <- list()
for(i in levels(DE.metadata$Cluster)){
  for(j in unique(DE.metadata$Cluster)){
    if(i != j){ # Don't compare to self
      if(! paste(j, "_over_", i, sep = "") %in% names(res.CEL.arch)){ # We don't need to do the reciprocals
        cat(paste(i, "_over_", j, "\n",sep = ""))
        res.CEL.arch[[paste(i, "_over_", j, sep = "")]] <- results(dds, contrast = c("Cluster", i, j), parallel = T)
        write.table(row.names(subset(res.CEL.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0)), file = paste("Pseudobulk.DE.CEL.RNA/archetype_results/", i, "_over_", j, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        write.table(row.names(subset(res.CEL.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0)), file = paste("Pseudobulk.DE.CEL.RNA/archetype_results/", j, "_over_", i, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        get_ontology_DE(res = subset(res.CEL.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0), name = paste(i, "_over_", j, sep = ""), outdir = "Pseudobulk.DE.CEL.RNA/archetype_results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
        get_ontology_DE(res = subset(res.CEL.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0), name = paste(j, "_over_", i, sep = ""), outdir = "Pseudobulk.DE.CEL.RNA/archetype_results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
      }
    }
  }
}

# Order the res tables on significance
for(i in names(res.CEL.arch)){
  res.CEL.arch[[i]] <- res.CEL.arch[[i]][order(res.CEL.arch[[i]]$padj),]
}

# Plot the top 50 DE genes
for(i in names(res.CEL.arch)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.CEL.arch[[i]])[1:50],])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.CEL.RNA/archetype_results/",i, ".top50_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.CEL.arch)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.CEL.arch[[i]], padj < 0.1)
  if(nrow(fc_ordered_res) < 1){next}
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors), 
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.CEL.RNA/archetype_results/",i, ".sig_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the top 50 DE genes - in the affected clusters only
for(i in names(res.CEL.arch)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.CEL.arch[[i]])[1:50],])
  a      <- strsplit(i, "_over_")[[1]][[1]]
  b      <- strsplit(i, "_over_")[[1]][[2]]
  keep   <- grep(a, colnames(de.mat))
  keep   <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.CEL.RNA/archetype_results/",i, ".top50_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.CEL.arch)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.CEL.arch[[i]], padj < 0.1)
  if(nrow(fc_ordered_res) < 1){next}
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  a              <- strsplit(i, "_over_")[[1]][[1]]
  b              <- strsplit(i, "_over_")[[1]][[2]]
  keep           <- grep(a, colnames(de.mat))
  keep           <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.CEL.RNA/archetype_results/",i, ".sig_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}


## Get DE genes for individual clusters -- Use 10X only RNA assay
dir.create("Pseudobulk.DE.10X.RNA/results", showWarnings = F, recursive = T)

# Remove monocyte clusters, look at the macs only
final.pop.call.from_full.integrated.mac.seurat <- subset(final.pop.call.from_full.integrated.mye.seurat, idents = as.vector(unique(Idents(final.pop.call.from_full.integrated.mye.seurat))[-grep("Monocytes",unique(Idents(final.pop.call.from_full.integrated.mye.seurat)))]))

# Remove CEL patients
final.pop.call.from_full.integrated.mac.10X_only.seurat <- subset(final.pop.call.from_full.integrated.mac.seurat, subset = Patient == "P1" | Patient == "P2" | Patient == "P3")
unique(final.pop.call.from_full.integrated.mac.10X_only.seurat$Patient)

# Prep a cluster + sample name grouping to make the pseudobulk counts over
cluster.patient <- paste(Idents(final.pop.call.from_full.integrated.mac.10X_only.seurat), final.pop.call.from_full.integrated.mac.10X_only.seurat$Patient, sep = ".")
cluster.patient <- gsub(pattern = " ", replacement = "_", x = cluster.patient)
cluster.patient <- gsub(pattern = "+", replacement = "hi_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "-", replacement = "lo_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Nonlo", replacement = "Non", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Monocytelo", replacement = "Monocyte_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "Residentlo", replacement = "Resident_", x = cluster.patient, fixed = T)
cluster.patient <- gsub(pattern = "__", replacement = "_", x = cluster.patient, fixed = T)

final.pop.call.from_full.integrated.mac.10X_only.seurat <- AddMetaData(final.pop.call.from_full.integrated.mac.10X_only.seurat, metadata = cluster.patient, col.name = "clus.sample")

# Make pseudobulk
final.pop.call.from_full.integrated.mac.avgexp <- AverageExpression(object = final.pop.call.from_full.integrated.mac.10X_only.seurat, group.by = "clus.sample", assays = "RNA")$RNA

# Make integers (to simulate readcounts)
final.pop.call.from_full.integrated.mac.avgexp <- ceiling(final.pop.call.from_full.integrated.mac.avgexp)

# Make metadata table
DE.metadata <- data.frame("Cluster" = t(as.data.frame(strsplit(colnames(final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,1], 
                          "Patient" = t(as.data.frame(strsplit(colnames(final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,2], 
                          row.names = colnames(final.pop.call.from_full.integrated.mac.avgexp))
colnames(DE.metadata) <- c("Cluster", "Patient")
DE.metadata$Cluster <- as.factor(DE.metadata$Cluster)
DE.metadata$Patient <- as.factor(DE.metadata$Patient)
head(DE.metadata)

# setup colors
pseudo.DE.colors <- M.int_refined.pop.colors
names(pseudo.DE.colors) <- gsub(pattern = " ", replacement = "_", x = names(pseudo.DE.colors))
names(pseudo.DE.colors) <- gsub(pattern = "+", replacement = "hi_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "-", replacement = "lo_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Nonlo", replacement = "Non", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Monocytelo", replacement = "Monocyte_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "Residentlo", replacement = "Resident_", x = names(pseudo.DE.colors), fixed = T)
names(pseudo.DE.colors) <- gsub(pattern = "__", replacement = "_", x = names(pseudo.DE.colors), fixed = T)
pseudo.DE.colors <- pseudo.DE.colors[-c(1,2,3)]


## Run DESeq2 for DE analyses
# Get the data ready
dds <- DESeqDataSetFromMatrix(countData = final.pop.call.from_full.integrated.mac.avgexp, colData = DE.metadata, design = ~ Patient + Cluster)

# Remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds, parallel = T)
vst <- varianceStabilizingTransformation(dds)
dds
vst

## Visualise general clustering
topVarGenes <- head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = DE.metadata,
         annotation_colors = list("Cluster" = pseudo.DE.colors),
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "Pseudobulk.DE.10X.RNA/top1000_variable_genes.kmeans_heatmap.pdf", width = 20, height = 20
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$Patient <- DE.metadata$Patient 
mat.pca$Cluster <- DE.metadata$Cluster
percentVar <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Cluster)) +
  geom_point(size=3) +
  scale_color_manual(values = pseudo.DE.colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.10X.RNA/top1000_variable_genes.PCA_plot_Cluster.pdf", width = 20, height = 10)

ggplot(mat.pca, aes(PC1, PC2, color = Patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.10X.RNA/top1000_variable_genes.PCA_plot_Patient.pdf")

resultsNames(dds)

res.10X <- list()
for(i in levels(DE.metadata$Cluster)){
  for(j in unique(DE.metadata$Cluster)){
    if(i != j){ # Don't compare to self
      if(! paste(j, "_over_", i, sep = "") %in% names(res.10X)){ # We don't need to do the reciprocals
        cat(paste(i, "_over_", j, "\n",sep = ""))
        res.10X[[paste(i, "_over_", j, sep = "")]] <- results(dds, contrast = c("Cluster", i, j), parallel = T)
        write.table(row.names(subset(res.10X[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0)), file = paste("Pseudobulk.DE.10X.RNA/results/", i, "_over_", j, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        write.table(row.names(subset(res.10X[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0)), file = paste("Pseudobulk.DE.10X.RNA/results/", j, "_over_", i, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        get_ontology_DE(res = subset(res.10X[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0), name = paste(i, "_over_", j, sep = ""), outdir = "Pseudobulk.DE.10X.RNA/results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
        get_ontology_DE(res = subset(res.10X[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0), name = paste(j, "_over_", i, sep = ""), outdir = "Pseudobulk.DE.10X.RNA/results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
      }
    }
  }
}

# Order the res tables on significance
for(i in names(res.10X)){
  res.10X[[i]] <- res.10X[[i]][order(res.10X[[i]]$padj),]
}

# Plot the top 50 DE genes
for(i in names(res.10X)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.10X[[i]])[1:50],])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.10X.RNA/results/",i, ".top50_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.10X)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.10X[[i]], padj < 0.1)
  if(nrow(fc_ordered_res) < 1){next}
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors), 
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.10X.RNA/results/",i, ".sig_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the top 50 DE genes - in the affected clusters only
for(i in names(res.10X)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.10X[[i]])[1:50],])
  a      <- strsplit(i, "_over_")[[1]][[1]]
  b      <- strsplit(i, "_over_")[[1]][[2]]
  keep   <- grep(a, colnames(de.mat))
  keep   <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.10X.RNA/results/",i, ".top50_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.10X)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.10X[[i]], padj < 0.1)
  if(nrow(fc_ordered_res) < 1){next}
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  a              <- strsplit(i, "_over_")[[1]][[1]]
  b              <- strsplit(i, "_over_")[[1]][[2]]
  keep           <- grep(a, colnames(de.mat))
  keep           <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.10X.RNA/results/",i, ".sig_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}


##===========================================================================
## Get DE genes for archetypes
dir.create("Pseudobulk.DE.10X.RNA/archetype_results", showWarnings = F, recursive = T)

## Set up archetypes
archetypes <- as.vector(Idents(final.pop.call.from_full.integrated.mye.seurat))
archetypes[grep("Foamy", archetypes)]        <- "Foamy"
archetypes[grep("Resident", archetypes)]     <- "Resident"
archetypes[grep("Lipid", archetypes)]        <- "Resident"
archetypes[grep("Monocyte-derived", archetypes)] <- "Inflammatory"

# Sanity check
levels(as.factor(archetypes))

# Add to the seurat object
final.pop.call.from_full.integrated.mye.seurat <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = archetypes, col.name = "archetype")

# Set up colors
archetype.colors <- c("Resident" = "#528347", "Inflammatory" = "#87482F", "Foamy" = "#9476A3", "Monocytes" = "#6495ED")

## Set up archetypes as default idents in a new object
archetype.integrated.mye.velocyto.seurat         <- AddMetaData(final.pop.call.from_full.integrated.mye.seurat, metadata = Idents(final.pop.call.from_full.integrated.mye.seurat), col.name = "final.pop.idents")
Idents(archetype.integrated.mye.velocyto.seurat) <- archetypes

# Remove monocyte clusters, look at the macs only
archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat <- subset(archetype.integrated.mye.velocyto.seurat, idents = as.vector(unique(Idents(archetype.integrated.mye.velocyto.seurat))[-grep("Monocytes",unique(Idents(archetype.integrated.mye.velocyto.seurat)))]))
unique(Idents(archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat))

# Remove CEL patients
archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat <- subset(archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat, subset = Patient == "P1" | Patient == "P2" | Patient == "P3")
unique(archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat$Patient)

# Prep a cluster + sample name grouping to make the pseudobulk counts over
cluster.patient <- paste(Idents(archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat), archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat$Patient, sep = ".")

archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat <- AddMetaData(archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat, metadata = cluster.patient, col.name = "clus.sample")

# Make pseudobulk
archetype.final.pop.call.from_full.integrated.mac.avgexp <- AverageExpression(object = archetype.final.pop.call.from_full.integrated.mac.10X_only.seurat, group.by = "clus.sample", assays = "RNA")$RNA

# Make integers (to simulate readcounts)
archetype.final.pop.call.from_full.integrated.mac.avgexp <- ceiling(archetype.final.pop.call.from_full.integrated.mac.avgexp)

# Make metadata table
DE.metadata <- data.frame("Cluster" = t(as.data.frame(strsplit(colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,1], 
                          "Patient" = t(as.data.frame(strsplit(colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp), ".", fixed = T)))[,2], 
                          row.names = colnames(archetype.final.pop.call.from_full.integrated.mac.avgexp))
colnames(DE.metadata) <- c("Cluster", "Patient")
DE.metadata$Cluster <- as.factor(DE.metadata$Cluster)
DE.metadata$Patient <- as.factor(DE.metadata$Patient)
head(DE.metadata)

# setup colors
pseudo.DE.colors <- archetype.colors

## Run DESeq2 for DE analyses
# Get the data ready
dds <- DESeqDataSetFromMatrix(countData = archetype.final.pop.call.from_full.integrated.mac.avgexp, colData = DE.metadata, design = ~ Patient + Cluster)

# Remove low count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

# Run DESeq
dds <- DESeq(dds, parallel = T)
vst <- varianceStabilizingTransformation(dds)
dds
vst

## Visualise general clustering
topVarGenes <- head(order(-rowVars(assay(vst))),1000) #Get 1000 most variable genes
mat <- assay(vst)[topVarGenes,]

#Kmeans heatmap
z.mat <- t(scale(t(mat)))
z.mat <- z.mat[is.finite(rowMeans(z.mat)),]
km.z.mat <- kmeans(z.mat,5,iter.max=20,nstart=20)
bk2 = unique(c(seq(min(z.mat), -0.01, length=63), 0, seq(0.01,max(z.mat), length=192)))
col1 = colorRampPalette(c("blue", "white"))(63)
col2 = colorRampPalette(c("white", "red"))(192)
mycols <- c(col1,"white",col2)

pheatmap(z.mat[names(sort(km.z.mat$cluster)),],
         cluster_rows = F,
         cluster_cols = T,
         breaks = bk2,
         annotation_row = as.data.frame(sort(km.z.mat$cluster)),
         annotation_col = DE.metadata,
         annotation_colors = list("Cluster" = pseudo.DE.colors),
         color = mycols,
         annotation_names_row = F,
         show_rownames = F,
         show_colnames = T,
         filename = "Pseudobulk.DE.10X.RNA/archetype_results/top1000_variable_genes.kmeans_heatmap.pdf", width = 20, height = 20
)

#2D PCA plot
mat.pca <- as.data.frame(prcomp(t(assay(vst)))$x)
mat.pca$Patient <- DE.metadata$Patient 
mat.pca$Cluster <- DE.metadata$Cluster
percentVar <- round(100 * prcomp(t(assay(vst)))$sdev^2/sum(prcomp(t(assay(vst)))$sdev^2))

ggplot(mat.pca, aes(PC1, PC2, color = Cluster)) +
  geom_point(size=3) +
  scale_color_manual(values = pseudo.DE.colors) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.10X.RNA/archetype_results/top1000_variable_genes.PCA_plot_Cluster.pdf", width = 20, height = 10)

ggplot(mat.pca, aes(PC1, PC2, color = Patient)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  ggtitle("PCA Plot") +
  theme_light() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.line = element_blank(),
        legend.position = "right"
  )
ggsave("Pseudobulk.DE.10X.RNA/archetype_results/top1000_variable_genes.PCA_plot_Patient.pdf")

resultsNames(dds)

res.10X.arch <- list()
for(i in levels(DE.metadata$Cluster)){
  for(j in unique(DE.metadata$Cluster)){
    if(i != j){ # Don't compare to self
      if(! paste(j, "_over_", i, sep = "") %in% names(res.10X.arch)){ # We don't need to do the reciprocals
        cat(paste(i, "_over_", j, "\n",sep = ""))
        res.10X.arch[[paste(i, "_over_", j, sep = "")]] <- results(dds, contrast = c("Cluster", i, j), parallel = T)
        write.table(row.names(subset(res.10X.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0)), file = paste("Pseudobulk.DE.10X.RNA/archetype_results/", i, "_over_", j, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        write.table(row.names(subset(res.10X.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0)), file = paste("Pseudobulk.DE.10X.RNA/archetype_results/", j, "_over_", i, ".sig.res.txt", sep = ""), quote = F, row.names = F, col.names = F)
        get_ontology_DE(res = subset(res.10X.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange > 0), name = paste(i, "_over_", j, sep = ""), outdir = "Pseudobulk.DE.10X.RNA/archetype_results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
        get_ontology_DE(res = subset(res.10X.arch[[paste(i, "_over_", j, sep = "")]], padj < 0.1 & log2FoldChange < 0), name = paste(j, "_over_", i, sep = ""), outdir = "Pseudobulk.DE.10X.RNA/archetype_results/", return.data = F, full_GSEA = F, plot.top.n = 10, volcano.plot = F, universe = final.pop.call.from_full.integrated.mye.seurat)
      }
    }
  }
}

# Order the res tables on significance
for(i in names(res.10X.arch)){
  res.10X.arch[[i]] <- res.10X.arch[[i]][order(res.10X.arch[[i]]$padj),]
}

# Plot the top 50 DE genes
for(i in names(res.10X.arch)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.10X.arch[[i]])[1:50],])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.10X.RNA/archetype_results/",i, ".top50_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.10X.arch)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.10X.arch[[i]], padj < 0.1)
  if(nrow(fc_ordered_res) < 1){next}
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  pheatmap(mat               = de.mat, 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors), 
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.10X.RNA/archetype_results/",i, ".sig_DE.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the top 50 DE genes - in the affected clusters only
for(i in names(res.10X.arch)){
  cat(paste(i, "\n",sep = ""))
  de.mat <- assay(vst[row.names(res.10X.arch[[i]])[1:50],])
  a      <- strsplit(i, "_over_")[[1]][[1]]
  b      <- strsplit(i, "_over_")[[1]][[2]]
  keep   <- grep(a, colnames(de.mat))
  keep   <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.10X.RNA/archetype_results/",i, ".top50_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

# Plot the sig DE genes
for(i in names(res.10X.arch)){
  cat(paste(i, "\n",sep = ""))
  fc_ordered_res <- subset(res.10X.arch[[i]], padj < 0.1)
  if(nrow(fc_ordered_res) < 1){next}
  fc_ordered_res <- fc_ordered_res[order(fc_ordered_res$log2FoldChange, decreasing = T),]
  de.mat         <- assay(vst[row.names(fc_ordered_res),])
  a              <- strsplit(i, "_over_")[[1]][[1]]
  b              <- strsplit(i, "_over_")[[1]][[2]]
  keep           <- grep(a, colnames(de.mat))
  keep           <- sort(unique(c(keep, grep(b, colnames(de.mat)))))
  pheatmap(mat               = de.mat[,keep], 
           scale             = "row",
           annotation_col    = DE.metadata,
           annotation_colors = list("Cluster" = pseudo.DE.colors),
           cluster_rows      = F,
           show_rownames     = T,
           show_colnames     = F,
           cellwidth         = 10, 
           cellheight        = 10, 
           drop_levels       = T,
           filename          = paste("Pseudobulk.DE.10X.RNA/archetype_results/",i, ".sig_DE.filtered.scaled_heatmap.pdf", sep = "")
  )
}

## Check the DE overlap bewteen the various methods and assays
res.arch.overlaps <- data.frame(row.names = names(res.10X.arch))
for(theIdent in names(res.10X.arch)){
  res.arch.overlaps[theIdent, 1] <- sum(row.names(subset(res.10X.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.CEL.arch[[theIdent]], padj < 0.1)))
  res.arch.overlaps[theIdent, 2] <- sum(row.names(subset(res.10X.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.RNA.arch[[theIdent]], padj < 0.1)))
  res.arch.overlaps[theIdent, 3] <- sum(row.names(subset(res.10X.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.arch.int[[theIdent]], padj < 0.1)))
  
  res.arch.overlaps[theIdent, 4] <- sum(row.names(subset(res.CEL.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.10X.arch[[theIdent]], padj < 0.1)))
  res.arch.overlaps[theIdent, 5] <- sum(row.names(subset(res.CEL.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.RNA.arch[[theIdent]], padj < 0.1)))
  res.arch.overlaps[theIdent, 6] <- sum(row.names(subset(res.CEL.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.arch.int[[theIdent]], padj < 0.1)))
  
  res.arch.overlaps[theIdent, 7] <- sum(row.names(subset(res.RNA.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.10X.arch[[theIdent]], padj < 0.1)))
  res.arch.overlaps[theIdent, 8] <- sum(row.names(subset(res.RNA.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.CEL.arch[[theIdent]], padj < 0.1)))
  res.arch.overlaps[theIdent, 9] <- sum(row.names(subset(res.RNA.arch[[theIdent]], padj < 0.1)) %in% row.names(subset(res.arch.int[[theIdent]], padj < 0.1)))
  
  res.arch.overlaps[theIdent, 10] <- sum(row.names(subset(res.arch.int[[theIdent]], padj < 0.1)) %in% row.names(subset(res.10X.arch[[theIdent]], padj < 0.1)))
  res.arch.overlaps[theIdent, 11] <- sum(row.names(subset(res.arch.int[[theIdent]], padj < 0.1)) %in% row.names(subset(res.RNA.arch[[theIdent]], padj < 0.1)))
  res.arch.overlaps[theIdent, 12] <- sum(row.names(subset(res.arch.int[[theIdent]], padj < 0.1)) %in% row.names(subset(res.CEL.arch[[theIdent]], padj < 0.1)))
}  

colnames(res.arch.overlaps) <- c("10X_in_CEL", "10X_in_RNA", "10X_in_INT", "CEL_in_10X", "CEL_in_RNA", "CEL_in_INT", "RNA_in_10X", "RNA_in_CEL", "RNA_in_INT", "INT_in_10X", "INT_in_RNA", "INT_in_CEL")
res.arch.overlaps

res.arch.totals <- data.frame(row.names = names(res.10X.arch))
for(theIdent in names(res.10X.arch)){
  res.arch.totals[theIdent, 1] <- nrow(subset(res.10X.arch[[theIdent]], padj < 0.1))
  res.arch.totals[theIdent, 2] <- nrow(subset(res.CEL.arch[[theIdent]], padj < 0.1))
  res.arch.totals[theIdent, 3] <- nrow(subset(res.arch.int[[theIdent]], padj < 0.1))
  res.arch.totals[theIdent, 4] <- nrow(subset(res.RNA.arch[[theIdent]], padj < 0.1))
}
colnames(res.arch.totals) <- c("10X", "CEL", "INT", "RNA")
res.arch.totals

res.arch.cors <- res.arch.overlaps[,1:3] / res.arch.totals$`10X`
res.arch.cors <- cbind(res.arch.cors, res.arch.overlaps[,4:6]   / res.arch.totals$`CEL`)
res.arch.cors <- cbind(res.arch.cors, res.arch.overlaps[,7:9]   / res.arch.totals$`RNA`)
res.arch.cors <- cbind(res.arch.cors, res.arch.overlaps[,10:12] / res.arch.totals$`INT`)
res.arch.cors

