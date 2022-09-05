## Export data for use with cellphonedb (pyhton)
dir.create("cellphonedb_input", showWarnings = F)
dir.create("cellphonedb_results", showWarnings = F)

# Subset out cells of interest (i.e., all plaque cells plus the PBMC monocytes)
LRI.cells <- WhichCells(integrated.full.seurat, expression = Tissue == "plaque")
LRI.cells <- c(LRI.cells, WhichCells(integrated.full.seurat, idents = levels(Idents(integrated.full.seurat))[grep("Monocytes", levels(Idents(integrated.full.seurat)))]))
LRI.cells.seurat <- subset(integrated.full.seurat, cells = LRI.cells)
customUMAP(LRI.cells.seurat, cols = full_set.colors, plot.width = 20, pt.size = 1, file.name = "cellphonedb_results/Cellphone input UMAP.pdf", legend.pos = "right")
customUMAP(LRI.cells.seurat, group.by = "Tissue", plot.width = 10, pt.size = 1, file.name = "cellphonedb_results/Cellphone input Tissue UMAP.pdf", legend.pos = "right")

# Counts
write.table(GetAssayData(LRI.cells.seurat, assay = "RNA", slot = "data"), file = 'cellphonedb_input/matrix.txt', quote = F , sep = "\t")

# Save gene names
write(x = rownames(GetAssayData(LRI.cells.seurat, assay = "RNA", slot = "data")), file = "cellphonedb_input/features.tsv")
write(x = colnames(GetAssayData(LRI.cells.seurat, assay = "RNA", slot = "data")), file = "cellphonedb_input/barcodes.tsv")

# Metadata
metaTable <- data.frame(Cell      = names(Idents(LRI.cells.seurat)), 
                        cell_type = Idents(LRI.cells.seurat))

unique(metaTable$cell_type)
head(metaTable)

# Remove spaces and '+' and '-' signs from the idents or downstream will be a mess. (cellphone will replace them with dots otherwise)
metaTable$cell_type <- gsub(x = metaTable$cell_type, pattern = " |\\+|-", replacement = "_")
metaTable$cell_type <- gsub(x = metaTable$cell_type, pattern = "__",    replacement = "_")
write.table(metaTable, file = "cellphonedb_input/metadata.table.txt", row.names = F, quote = F, sep = "\t")

# Make a conversion table for idents
cell_type_map <- data.frame(cpdb.ident = unique(metaTable$cell_type), seurat.ident = unique(Idents(LRI.cells.seurat)))

# Marker genes
# Run FindAllMarkers again to update it with the refined clusters, if needed
LRI.cells.seurat.markers <- FindAllMarkers(LRI.cells.seurat, assay = "RNA", only.pos = TRUE, min.pct = 0.25, thresh = 0.25)

# Keep sig DEGs, and clean up the names again
fDEGs         <- subset(LRI.cells.seurat.markers, p_val_adj < 0.1)
fDEGs$cluster <- gsub(x = fDEGs$cluster, pattern = " |\\+|-|\\.", replacement = "_")
fDEGs$cluster <- gsub(x = fDEGs$cluster, pattern = "__",    replacement = "_")

# 1st column = cluster; 2nd column = gene 
fDEGs <- fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
write.table(fDEGs, file ='cellphonedb_input/cellphonedb.marker.input.txt', sep = '\t', quote = F, row.names = F)
head(fDEGs)
unique(fDEGs$cluster)

#===========================================================================================================================
#===========================================================================================================================
dir.create("cellphonedb_results/lig.rec.int", showWarnings = F)

# Import results from cellphonedb
cpdb.means <- read.table(file = "cellphonedb_input/DEG_analysis_out/significant_means.txt", row.names = 1, header = T, sep = "\t")

# Keep only relevant interactions (involving at least 1 DEG)
cpdb.relevant <- read.table(file = "cellphonedb_input/DEG_analysis_out/relevant_interactions.txt", row.names = 1, header = T, sep = "\t")
cpdb.relevant <- cpdb.means[row.names(cpdb.means) %in% row.names(cpdb.relevant),]

# Get deconvolution of complexes
cpdb.deconv <- read.table(file = "cellphonedb_input/DEG_analysis_out/deconvoluted.txt", header = T, sep = "\t")
cpdb.deconv <- cpdb.deconv[cpdb.deconv$id_cp_interaction %in% row.names(cpdb.relevant),]

cpdb.relevant[is.na(cpdb.relevant)] <- 0
cpdb.relevant[1:5,1:15]

# Remove integrin interactions (meh let's keep 'em for now)
# cpdb.relevant <- cpdb.relevant[cpdb.relevant$is_integrin == "False",]

#===========================================================================================================================
# Plot overall interactions
# Shape the data for plotting
s.cpdb.relevant <- cpdb.relevant
s.cpdb.relevant[s.cpdb.relevant$receptor_a == "True", "receptor_a"] <- TRUE
s.cpdb.relevant[s.cpdb.relevant$receptor_a == "False","receptor_a"] <- FALSE
s.cpdb.relevant[s.cpdb.relevant$receptor_b == "True", "receptor_b"] <- TRUE
s.cpdb.relevant[s.cpdb.relevant$receptor_b == "False","receptor_b"] <- FALSE
s.cpdb.relevant[1:10,1:12]

# Prep ligand and receptor name columns
s.cpdb.relevant$ligand   <- "NA"
s.cpdb.relevant$receptor <- "NA"

#Eek out ligands and receptors using the 'is receptor' and gen info
attach(s.cpdb.relevant)

# First deal with the single gene entries
# gene A is a receptor and not a complex
index <- which(receptor_a == TRUE & gene_a != "")
s.cpdb.relevant[index, "receptor"] <- s.cpdb.relevant[index,"gene_a"]

# genes A is a receptor and its ligand gene b is not a complex
index <- which(receptor_a == TRUE & gene_b != "")
s.cpdb.relevant[index, "ligand"] <- s.cpdb.relevant[index,"gene_b"]

# gene B is a receptor and not a complex
index <- which(receptor_b == TRUE & gene_b != "")
s.cpdb.relevant[index, "receptor"] <- s.cpdb.relevant[index,"gene_b"]

# genes B is a receptor and its ligand gene a is not a complex
index <- which(receptor_b == TRUE & gene_a != "")
s.cpdb.relevant[index, "ligand"] <- s.cpdb.relevant[index,"gene_a"]

# Now do the complex entries
# gene A is a receptor and a complex
index <- which(receptor_a == TRUE & gene_a == "")
s.cpdb.relevant[index, "receptor"] <- s.cpdb.relevant[index,"partner_a"]

# genes A is a receptor and its ligand gene b is a complex
index <- which(receptor_a == TRUE & gene_b == "")
s.cpdb.relevant[index, "ligand"] <- s.cpdb.relevant[index,"partner_b"]

# gene B is a receptor and a complex
index <- which(receptor_b == TRUE & gene_b == "")
s.cpdb.relevant[index, "receptor"] <- s.cpdb.relevant[index,"partner_b"]

# genes B is a receptor and its ligand gene a is a complex
index <- which(receptor_b == TRUE & gene_a == "")
s.cpdb.relevant[index, "ligand"] <- s.cpdb.relevant[index,"partner_a"]

# Remove the 'no recptor pairs'
index <-  which(receptor_a == FALSE & receptor_b == FALSE)
s.cpdb.relevant <- s.cpdb.relevant[-index,]
detach(s.cpdb.relevant)

# Melt the 'frame using ligand and receptor info and the interaction scores (take -2 the length of the colname sbecause we appendend the ligand adn receptor columns to the end of the df)
m.s.cpdb.relevant <- melt(s.cpdb.relevant, id.vars = c("ligand", "receptor"), measure.vars = c(12:(length(colnames(s.cpdb.relevant))-2)), variable.name = "cells", value.name = "interaction.score")
tail(m.s.cpdb.relevant)

# Subset pairs with significant interactions (interaction score > 0)
exp.m.s.cpdb.relevant <- m.s.cpdb.relevant[m.s.cpdb.relevant$interaction.score > 0,]
dim(exp.m.s.cpdb.relevant)

# Prep ligand and receptor cell names columns
l <- vector()
r <- vector()
for (i in strsplit(as.character(exp.m.s.cpdb.relevant$cells), split = ".", fixed = T)){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
exp.m.s.cpdb.relevant$ligand.cells <- l
exp.m.s.cpdb.relevant$receptor.cells <- r

nrow(exp.m.s.cpdb.relevant)
head(exp.m.s.cpdb.relevant)
tail(exp.m.s.cpdb.relevant)

# Clean up
exp.m.s.cpdb.relevant$interaction.score <- as.numeric(exp.m.s.cpdb.relevant$interaction.score)
exp.m.s.cpdb.relevant <- exp.m.s.cpdb.relevant[!is.na(exp.m.s.cpdb.relevant$interaction.score),]

# Fix factors
exp.m.s.cpdb.relevant$receptor       <- droplevels(as.factor(exp.m.s.cpdb.relevant$receptor))
exp.m.s.cpdb.relevant$ligand         <- droplevels(as.factor(exp.m.s.cpdb.relevant$ligand))
exp.m.s.cpdb.relevant$receptor.cells <- droplevels(as.factor(exp.m.s.cpdb.relevant$receptor.cells))
exp.m.s.cpdb.relevant$ligand.cells   <- droplevels(as.factor(exp.m.s.cpdb.relevant$ligand.cells))
length(levels(exp.m.s.cpdb.relevant$receptor.cells))
levels(exp.m.s.cpdb.relevant$receptor.cells)

# Define colors for the subclusters
subclus.colors <- full_set.colors

# Sort approximately the same as the levels in our df for convenience, then order manually
subclus.colors <- subclus.colors[sort(names(subclus.colors))]
subclus.colors <- subclus.colors[names(subclus.colors)[c(1,3,5,4,7,2,8,9,6,10,14,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28)]]
names(subclus.colors) <- levels(exp.m.s.cpdb.relevant$ligand.cells)

# Plot all interactions
ggplot(exp.m.s.cpdb.relevant, aes(x = ligand, y= interaction.score, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.ligands_in_receptor_cells.outliers_removed.subset.dots.pdf", width = 40, height = 10)

ggplot(exp.m.s.cpdb.relevant, aes(x = receptor, y= interaction.score, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.receptors_in_ligand_cells.outliers_removed.subset.dots.pdf", width = 40, height = 10)

ggplot(exp.m.s.cpdb.relevant, aes(x = ligand, y= interaction.score, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.ligands_in_ligand_cells.outliers_removed.subset.dots.pdf", width = 40, height = 10)

ggplot(exp.m.s.cpdb.relevant, aes(x = receptor, y= interaction.score, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.receptors_in_receptor_cells.outliers_removed.subset.dots.pdf", width = 40, height = 10)

#=============================================================================================================
# PLot only from specific population (group) viewpoints
#=============================================================================================================
## First just do from all to all populations
dir.create( "cellphonedb_results/per_population", showWarnings = F, recursive = T)
for(theIdent in unique(c(exp.m.s.cpdb.relevant$ligand.cells, exp.m.s.cpdb.relevant$receptor.cells))){
  cat(paste("Working on: ", theIdent, "...\n", sep = ""))
  custom_LRI_plots(save.dir = "cellphonedb_results/per_population", name = theIdent, sub.re = theIdent, subclus.colors = subclus.colors, verbose = T, top.score = 1)
}

# Now check the myeloid populations
custom_LRI_plots(save.dir = "cellphonedb_results/", name = "Myeloid", sub.re = "Mono|Macro", subclus.colors = subclus.colors, verbose = T, top.score = 1)

# Now check the resident populations
custom_LRI_plots(save.dir = "cellphonedb_results/", name = "Resident", sub.re = "TREM2_FOLR2", subclus.colors = subclus.colors, verbose = T, top.score = 1)

# Now check the Foamy populations
custom_LRI_plots(save.dir = "cellphonedb_results/", name = "Foamy", sub.re = "Foamy", subclus.colors = subclus.colors, verbose = T, top.score = 1)

# Now check the Inflammatory populations
custom_LRI_plots(save.dir = "cellphonedb_results/", name = "Inflammatory", sub.re = "IL1B_SELL", subclus.colors = subclus.colors, verbose = T, top.score = 1)

# Now check the monocyte populations
custom_LRI_plots(save.dir = "cellphonedb_results/", name = "Monocytic", sub.re = "Monocytes", subclus.colors = subclus.colors, verbose = T, top.score = 1)


## Let's extract the lists per group also
tail(exp.m.s.cpdb.relevant)
## Per Population
# Ligands
for(theIdent in unique(exp.m.s.cpdb.relevant$ligand.cells)){
  write.table(subset(exp.m.s.cpdb.relevant, ligand.cells == theIdent), file = paste("cellphonedb_results/per_population/", theIdent, "_ligand.interactions.txt", sep  = ""), 
              row.names = F, quote = F, sep = "\t")
}
# Receptors
for(theIdent in unique(exp.m.s.cpdb.relevant$receptor.cells)){
  write.table(subset(exp.m.s.cpdb.relevant, receptor.cells == theIdent), file = paste("cellphonedb_results/per_population/", theIdent, "_receptor.interactions.txt", sep  = ""),
              row.names = F, quote = F, sep = "\t")
}

## Myeloid populations
# Ligands
write.table(exp.m.s.cpdb.relevant[grep("Mono|Macro", exp.m.s.cpdb.relevant$ligand.cells),], file = "cellphonedb_results/Myeloid_cells_ligand.interactions.txt", 
            row.names = F, quote = F, sep = "\t")
# Receptors
write.table(exp.m.s.cpdb.relevant[grep("Mono|Macro", exp.m.s.cpdb.relevant$receptor.cells),], file = "cellphonedb_results/Myeloid_cells_receptor.interactions.txt", 
            row.names = F, quote = F, sep = "\t")

## Resident populations
# Ligands
write.table(exp.m.s.cpdb.relevant[grep("TREM2_FOLR", exp.m.s.cpdb.relevant$ligand.cells),], file = "cellphonedb_results/Resident_cells_ligand.interactions.txt", 
            row.names = F, quote = F, sep = "\t")
# Receptors
write.table(exp.m.s.cpdb.relevant[grep("TREM2_FOLR", exp.m.s.cpdb.relevant$receptor.cells),], file = "cellphonedb_results/Resident_cells_receptor.interactions.txt", 
            row.names = F, quote = F, sep = "\t")

## Foamy populations
# Ligands
write.table(exp.m.s.cpdb.relevant[grep("Foamy", exp.m.s.cpdb.relevant$ligand.cells),], file = "cellphonedb_results/Foamy_cells_ligand.interactions.txt", 
            row.names = F, quote = F, sep = "\t")
# Receptors
write.table(exp.m.s.cpdb.relevant[grep("Foamy", exp.m.s.cpdb.relevant$receptor.cells),], file = "cellphonedb_results/Foamy_cells_receptor.interactions.txt", 
            row.names = F, quote = F, sep = "\t")

## Inflammatory populations
# Ligands
write.table(exp.m.s.cpdb.relevant[grep("IL1B_SELL", exp.m.s.cpdb.relevant$ligand.cells),], file = "cellphonedb_results/Inflammatory_cells_ligand.interactions.txt", 
            row.names = F, quote = F, sep = "\t")
# Receptors
write.table(exp.m.s.cpdb.relevant[grep("IL1B_SELL", exp.m.s.cpdb.relevant$receptor.cells),], file = "cellphonedb_results/Inflammatory_cells_receptor.interactions.txt", 
            row.names = F, quote = F, sep = "\t")

## Monocyte populations
# Ligands
write.table(exp.m.s.cpdb.relevant[grep("Monocytes", exp.m.s.cpdb.relevant$ligand.cells),], file = "cellphonedb_results/Monocytic_cells_ligand.interactions.txt", 
            row.names = F, quote = F, sep = "\t")
# Receptors
write.table(exp.m.s.cpdb.relevant[grep("Monocytes", exp.m.s.cpdb.relevant$receptor.cells),], file = "cellphonedb_results/Monocytic_cells_receptor.interactions.txt", 
            row.names = F, quote = F, sep = "\t")



## Plot a heatmap of the number of interactions per celltype pair
# Get the number of interactions per cell type pair
count.cpdb.relevant <- cpdb.relevant[,12:ncol(cpdb.relevant)] > 0
count.cpdb.relevant <- colSums(count.cpdb.relevant)

## Shape the data for plotting
# Split the cell pairs
m.count.cpdb.relevant            <- melt(count.cpdb.relevant)
m.count.cpdb.relevant$A          <- sub("\\..*", "", row.names(m.count.cpdb.relevant))
m.count.cpdb.relevant$B          <- sub(".*.\\.", "", row.names(m.count.cpdb.relevant))
row.names(m.count.cpdb.relevant) <- NULL

# Turn into a df
mat.count.cpdb.relevant <- data.frame(row.names =  unique(m.count.cpdb.relevant$A))
for(celltype in unique(m.count.cpdb.relevant$A)){
  mat.count.cpdb.relevant[,celltype] <- m.count.cpdb.relevant[m.count.cpdb.relevant$A == celltype,"value"]
}

# Sanity check
mat.count.cpdb.relevant[1:5,1:5]
head(m.count.cpdb.relevant)

# Order on column of highest count
high.count.cells        <- colnames(mat.count.cpdb.relevant)[apply(mat.count.cpdb.relevant, 1, max) == max(mat.count.cpdb.relevant)]
mat.count.cpdb.relevant <- mat.count.cpdb.relevant[order(mat.count.cpdb.relevant[,high.count.cells], decreasing = T),]
mat.count.cpdb.relevant <- mat.count.cpdb.relevant[,row.names(mat.count.cpdb.relevant)]

## Make the plot
pheatmap(mat = mat.count.cpdb.relevant, color = colorSpacer(startcolor = "black", endcolor = "red", steps = max(mat.count.cpdb.relevant)), border_color = NA,
         cellwidth = 20, cellheight = 20, filename = "cellphonedb_results/heatmap.pdf", cluster_rows = F, cluster_cols = F, width = 15, height = 15)
