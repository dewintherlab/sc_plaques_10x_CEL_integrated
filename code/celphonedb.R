## Export data for use with cellphonedb (pyhton)
dir.create("cellphonedb_input", showWarnings = F)

# Counts
write.table(GetAssayData(integrated.full.seurat, assay = "RNA", slot = "data"), file = 'cellphonedb_input/matrix.txt', quote = F , sep = "\t")

# Save gene names
write(x = rownames(GetAssayData(integrated.full.seurat, assay = "RNA", slot = "data")), file = "cellphonedb_input/features.tsv")
write(x = colnames(GetAssayData(integrated.full.seurat, assay = "RNA", slot = "data")), file = "cellphonedb_input/barcodes.tsv")

# Metadata
metaTable <- data.frame(Cell      = names(Idents(integrated.full.seurat)), 
                        cell_type = Idents(integrated.full.seurat))

unique(metaTable$cell_type)
head(metaTable)

# Remove spaces and '+' and '-' signs from the idents or downstream will be a mess. (cellphone will replace them with dots otherwise)
metaTable$cell_type <- gsub(x = metaTable$cell_type, pattern = " |\\+|-", replacement = "_")
metaTable$cell_type <- gsub(x = metaTable$cell_type, pattern = "__",    replacement = "_")
write.table(metaTable, file = "cellphonedb_input/metadata.table.txt", row.names = F, quote = F, sep = "\t")

# Marker genes
# Run FindAllMarkers again to update it with the refined clusters, if needed
integrated.full.seurat.markers <- FindAllMarkers(integrated.full.seurat, assay = "RNA", only.pos = TRUE, min.pct = 0.25, thresh = 0.25)

# Keep sig DEGs, and clean up the names again
fDEGs         <- subset(integrated.full.seurat.markers, p_val_adj < 0.1)
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
cpdb.means <- read.table(file = "cellphonedb_results/significant_means.txt", row.names = 1, header = T, sep = "\t")

# Keep only relevant interactions (involving at least 1 DEG)
cpdb.relevant <- read.table(file = "cellphonedb_results/relevant_interactions.txt", row.names = 1, header = T, sep = "\t")
cpdb.relevant <- cpdb.means[row.names(cpdb.means) %in% row.names(cpdb.relevant),]

# Get deconvolution of complexes
cpdb.deconv <- read.table(file = "cellphonedb_results/deconvoluted.txt", header = T, sep = "\t")
cpdb.deconv <- cpdb.deconv[cpdb.deconv$id_cp_interaction %in% row.names(cpdb.relevant),]

cpdb.relevant[is.na(cpdb.relevant)] <- 0
cpdb.relevant[1:5,1:15]

# Remove integrin interactions
cpdb.relevant <- cpdb.relevant[cpdb.relevant$is_integrin == "False",]

# Save all HLA interactions for later
hla.cpdb.relevant <-cpdb.relevant[grep("HLA",cpdb.relevant$interacting_pair),]

# Keep only top ranked interactions
summary(cpdb.relevant$rank)
cpdb.relevant <- cpdb.relevant[cpdb.relevant$rank < mean(cpdb.relevant$rank),]
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
detach(s.cpdb.relevant)


# Melt the 'frame using ligand and receptor info and the interaction scores (take -2 the length of the colname sbecause we appendend the ligand adn receptor columns to the end of the df)
m.s.cpdb.relevant <- melt(s.cpdb.relevant, id.vars = c("ligand", "receptor"), measure.vars = c(12:(length(colnames(s.cpdb.relevant))-2)), variable.name = "cells", value.name = "interaction.score")
tail(m.s.cpdb.relevant)

# Subset pairs with significant interactions (>0)
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

#Define colors for the subclusters
subclus.colors <- c( "antiquewhite3", # ACTA2
                     colorSpacer(startcolor = "goldenrod1", endcolor = "goldenrod4", steps = 6, return.colors = T), # CD3
                     "darkolivegreen3", "darkolivegreen4", # Endo
                     "darkorange1", "darkorange3", # NK
                     colorSpacer(startcolor = "deepskyblue", endcolor = "blue4", steps = 7), # fullloid
                     "deeppink", # mast
                     "darkorchid1", "darkorchid3", # B
                     "coral4") # FOXP3
names(subclus.colors) <- levels(exp.m.s.cpdb.relevant$ligand.cells)

# Subset a bit for clarity
summary(exp.m.s.cpdb.relevant$interaction.score)
sub.exp.m.s.cpdb.relevant <- exp.m.s.cpdb.relevant[exp.m.s.cpdb.relevant$interaction.score > summary(exp.m.s.cpdb.relevant$interaction.score)[5],]
sub.exp.m.s.cpdb.relevant <- exp.m.s.cpdb.relevant[exp.m.s.cpdb.relevant$interaction.score < 20,]

# plot
ggplot(sub.exp.m.s.cpdb.relevant, aes(x = ligand, y= interaction.score, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.ligands_in_receptor_cells.outliers_removed.subset.dots.pdf", width = 30, height = 10)

ggplot(sub.exp.m.s.cpdb.relevant, aes(x = receptor, y= interaction.score, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.receptors_in_ligand_cells.outliers_removed.subset.dots.pdf", width = 30, height = 10)

ggplot(sub.exp.m.s.cpdb.relevant, aes(x = ligand, y= interaction.score, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.ligands_in_ligand_cells.outliers_removed.subset.dots.pdf", width = 30, height = 10)

ggplot(sub.exp.m.s.cpdb.relevant, aes(x = receptor, y= interaction.score, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.receptors_in_receptor_cells.outliers_removed.subset.dots.pdf", width = 30, height = 10)

#===========================================================================================================================
# Subset fullloid ligand interactions
full.lig.cells <- levels(exp.m.s.cpdb.relevant$ligand.cells)[grep("mac|Foam", levels(exp.m.s.cpdb.relevant$ligand.cells))]
exp.m.sig.full.lig.cpdb.relevant <- exp.m.s.cpdb.relevant[exp.m.s.cpdb.relevant$ligand.cells %in% full.lig.cells, ]
head(exp.m.sig.full.lig.cpdb.relevant)

# Check dimensions
length(unique(exp.m.sig.full.lig.cpdb.relevant$ligand))
length(unique(exp.m.sig.full.lig.cpdb.relevant$receptor))

# And plot
ggplot(exp.m.sig.full.lig.cpdb.relevant, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("fullloid ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in fullloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("cellphonedb_results/fullloid.ligands.dots.pdf", width = 25, height = 25)

# Subset a bit more for clarity
summary(exp.m.sig.full.lig.cpdb.relevant$interaction.score)
sub.exp.m.sig.full.lig.cpdb.relevant <- exp.m.sig.full.lig.cpdb.relevant[exp.m.sig.full.lig.cpdb.relevant$interaction.score > summary(exp.m.sig.full.lig.cpdb.relevant$interaction.score)[5],]

summary(summary(exp.m.sig.full.lig.cpdb.relevant$ligand))
sub.exp.m.sig.full.lig.cpdb.relevant <- sub.exp.m.sig.full.lig.cpdb.relevant[sub.exp.m.sig.full.lig.cpdb.relevant$ligand %in% names(summary(sub.exp.m.sig.full.lig.cpdb.relevant$ligand))[summary(summary(sub.exp.m.sig.full.lig.cpdb.relevant$ligand)) < 10],]

# Check dimensions
length(unique(sub.exp.m.sig.full.lig.cpdb.relevant$ligand))
length(unique(sub.exp.m.sig.full.lig.cpdb.relevant$receptor))

# And plot
ggplot(sub.exp.m.sig.full.lig.cpdb.relevant, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("fullloid ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in fullloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 14),
        aspect.ratio = 1)
ggsave("cellphonedb_results/fullloid.ligands.subset.dots.pdf", width = 15, height = 15)


#===========================================================================================================================
# Subset fullloid receptor interactions
full.rec.cells <- levels(exp.m.s.cpdb.relevant$receptor.cells)[grep("mac|Foam", levels(exp.m.s.cpdb.relevant$receptor.cells))]
exp.m.sig.full.rec.cpdb.relevant <- exp.m.s.cpdb.relevant[exp.m.s.cpdb.relevant$receptor.cells %in% full.rec.cells, ]
head(exp.m.sig.full.rec.cpdb.relevant)

# Check dimensions
length(unique(exp.m.sig.full.rec.cpdb.relevant$ligand))
length(unique(exp.m.sig.full.rec.cpdb.relevant$receptor))

# And plot
ggplot(exp.m.sig.full.rec.cpdb.relevant, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("fullloid receptors")) +
  theme_light() +
  ylab("Receptors (in fullloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("cellphonedb_results/fullloid.receptors.dots.pdf", width = 25, height = 25)

# Subset a bit more for clarity
summary(exp.m.sig.full.rec.cpdb.relevant$interaction.score)
sub.exp.m.sig.full.rec.cpdb.relevant <- exp.m.sig.full.rec.cpdb.relevant[exp.m.sig.full.rec.cpdb.relevant$interaction.score > summary(exp.m.sig.full.rec.cpdb.relevant$interaction.score)[5],]

summary(summary(exp.m.sig.full.rec.cpdb.relevant$ligand))
sub.exp.m.sig.full.rec.cpdb.relevant <- sub.exp.m.sig.full.rec.cpdb.relevant[sub.exp.m.sig.full.rec.cpdb.relevant$ligand %in% names(summary(sub.exp.m.sig.full.rec.cpdb.relevant$ligand))[summary(sub.exp.m.sig.full.rec.cpdb.relevant$ligand) < 20],]

# Check dimensions
length(unique(sub.exp.m.sig.full.rec.cpdb.relevant$ligand))
length(unique(sub.exp.m.sig.full.rec.cpdb.relevant$receptor))

# And plot
ggplot(sub.exp.m.sig.full.rec.cpdb.relevant, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("fullloid receptors")) +
  theme_light() +
  ylab("Receptors (in fullloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("cellphonedb_results/fullloid.receptors.subset.dots.pdf", width = 15, height = 15)




#=============================================================================================================
# Investigate random interactions

sub.exp.m.s.cpdb.relevant[grep("FAM3C", sub.exp.m.s.cpdb.relevant$receptor),][sub.exp.m.s.cpdb.relevant[grep("FAM3C", sub.exp.m.s.cpdb.relevant$receptor),"interaction.score"] > 15,]


#=============================================================================================================
# PLot HLA interactions in general and in and with fullloid cells
#=============================================================================================================
#===========================================================================================================================
# Plot overall interactions
# Shape the data for plotting
s.hla.spdb.means <- hla.cpdb.relevant
s.hla.spdb.means[s.hla.spdb.means$receptor_a == "True", "receptor_a"] <- TRUE
s.hla.spdb.means[s.hla.spdb.means$receptor_a == "False","receptor_a"] <- FALSE
s.hla.spdb.means[s.hla.spdb.means$receptor_b == "True", "receptor_b"] <- TRUE
s.hla.spdb.means[s.hla.spdb.means$receptor_b == "False","receptor_b"] <- FALSE
s.hla.spdb.means[1:10,1:12]

# Prep ligand and receptor name columns
s.hla.spdb.means$ligand   <- "NA"
s.hla.spdb.means$receptor <- "NA"

#Eek out ligands and receptors using the 'is receptor' and gen info
attach(s.hla.spdb.means)

# First deal with the single gene entries
# gene A is a receptor and not a complex
index <- which(receptor_a == TRUE & gene_a != "")
s.hla.spdb.means[index, "receptor"] <- s.hla.spdb.means[index,"gene_a"]

# genes A is a receptor and its ligand gene b is not a complex
index <- which(receptor_a == TRUE & gene_b != "")
s.hla.spdb.means[index, "ligand"] <- s.hla.spdb.means[index,"gene_b"]

# gene B is a receptor and not a complex
index <- which(receptor_b == TRUE & gene_b != "")
s.hla.spdb.means[index, "receptor"] <- s.hla.spdb.means[index,"gene_b"]

# genes B is a receptor and its ligand gene a is not a complex
index <- which(receptor_b == TRUE & gene_a != "")
s.hla.spdb.means[index, "ligand"] <- s.hla.spdb.means[index,"gene_a"]

# Now do the complex entries
# gene A is a receptor and a complex
index <- which(receptor_a == TRUE & gene_a == "")
s.hla.spdb.means[index, "receptor"] <- s.hla.spdb.means[index,"partner_a"]

# genes A is a receptor and its ligand gene b is a complex
index <- which(receptor_a == TRUE & gene_b == "")
s.hla.spdb.means[index, "ligand"] <- s.hla.spdb.means[index,"partner_b"]

# gene B is a receptor and a complex
index <- which(receptor_b == TRUE & gene_b == "")
s.hla.spdb.means[index, "receptor"] <- s.hla.spdb.means[index,"partner_b"]

# genes B is a receptor and its ligand gene a is a complex
index <- which(receptor_b == TRUE & gene_a == "")
s.hla.spdb.means[index, "ligand"] <- s.hla.spdb.means[index,"partner_a"]
detach(s.hla.spdb.means)


# Melt the 'frame using ligand and receptor info and the interaction scores (take -2 the length of the colname sbecause we appendend the ligand adn receptor columns to the end of the df)
m.s.hla.spdb.means <- melt(s.hla.spdb.means, id.vars = c("ligand", "receptor"), measure.vars = c(12:(length(colnames(s.hla.spdb.means))-2)), variable.name = "cells", value.name = "interaction.score")
tail(m.s.hla.spdb.means)

# Subset pairs with significant interactions (>0)
exp.m.s.hla.spdb.means <- m.s.hla.spdb.means[m.s.hla.spdb.means$interaction.score > 0,]
dim(exp.m.s.hla.spdb.means)

# Prep ligand and receptor cell names columns
l <- vector()
r <- vector()
exp.m.s.hla.spdb.means$cells <- gsub(pattern = "ass.swi", replacement = "ass_swi", x = exp.m.s.hla.spdb.means$cells, fixed = T)
for (i in strsplit(as.character(exp.m.s.hla.spdb.means$cells), split = ".", fixed = T)){
  l <- c(l, i[1])
  r <- c(r, i[2])
}
exp.m.s.hla.spdb.means$ligand.cells <- l
exp.m.s.hla.spdb.means$receptor.cells <- r

nrow(exp.m.s.hla.spdb.means)
head(exp.m.s.hla.spdb.means)
tail(exp.m.s.hla.spdb.means)

# Clean up
exp.m.s.hla.spdb.means$interaction.score <- as.numeric(exp.m.s.hla.spdb.means$interaction.score)
exp.m.s.hla.spdb.means <- exp.m.s.hla.spdb.means[!is.na(exp.m.s.hla.spdb.means$interaction.score),]

# Remove 'mixed cells' and 'My.4' sub clusters from the receptor and ligand cells as they are uninformative.
#exp.m.s.hla.spdb.means <- exp.m.s.hla.spdb.means[grep("My.4|Mixed", exp.m.s.hla.spdb.means$receptor.cells, invert = T),]
#exp.m.s.hla.spdb.means <- exp.m.s.hla.spdb.means[grep("My.4|Mixed", exp.m.s.hla.spdb.means$ligand.cells, invert = T),]

# Fix factors
exp.m.s.hla.spdb.means$receptor       <- droplevels(as.factor(exp.m.s.hla.spdb.means$receptor))
exp.m.s.hla.spdb.means$ligand         <- droplevels(as.factor(exp.m.s.hla.spdb.means$ligand))
exp.m.s.hla.spdb.means$receptor.cells <- droplevels(as.factor(exp.m.s.hla.spdb.means$receptor.cells))
exp.m.s.hla.spdb.means$ligand.cells   <- droplevels(as.factor(exp.m.s.hla.spdb.means$ligand.cells))
length(levels(exp.m.s.hla.spdb.means$receptor.cells))
levels(exp.m.s.hla.spdb.means$receptor.cells)

#Define colors for the subclusters
subclus.colors <- c( "antiquewhite3", # ACTA2
                     colorSpacer(startcolor = "goldenrod1", endcolor = "goldenrod4", steps = 6, return.colors = T), # CD3
                     "darkolivegreen3", "darkolivegreen4", # Endo
                     "darkorange1", "darkorange3", # NK
                     colorSpacer(startcolor = "deepskyblue", endcolor = "blue4", steps = 5), # fullloid
                     "deeppink", # mast
                     "darkorchid1", "darkorchid3", # B
                     "coral4") # FOXP3
names(subclus.colors) <- levels(exp.m.s.hla.spdb.means$ligand.cells)

# Subset a bit for clarity
summary(exp.m.s.hla.spdb.means$interaction.score)
sub.exp.m.s.hla.spdb.means <- exp.m.s.hla.spdb.means[exp.m.s.hla.spdb.means$interaction.score > summary(exp.m.s.hla.spdb.means$interaction.score)[5],]
sub.exp.m.s.hla.spdb.means <- exp.m.s.hla.spdb.means[exp.m.s.hla.spdb.means$interaction.score < 20,]

# plot
ggplot(sub.exp.m.s.hla.spdb.means, aes(x = ligand, y= interaction.score, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/HLA.All.ligands_in_receptor_cells.outliers_removed.subset.dots.pdf", width = 10, height = 10)

ggplot(sub.exp.m.s.hla.spdb.means, aes(x = receptor, y= interaction.score, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/HLA>All.receptors_in_ligand_cells.outliers_removed.subset.dots.pdf", width = 10, height = 10)

ggplot(sub.exp.m.s.hla.spdb.means, aes(x = ligand, y= interaction.score, col = ligand.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("ligands in ligand cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.ligands_in_ligand_cells.outliers_removed.subset.dots.pdf", width = 10, height = 10)

ggplot(sub.exp.m.s.hla.spdb.means, aes(x = receptor, y= interaction.score, col = receptor.cells)) +
  geom_point(size = 2) +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("receptors in receptor cells")) +
  ylab("Interaction score") +
  theme_light() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1/5
  )
ggsave("cellphonedb_results/All.receptors_in_receptor_cells.outliers_removed.subset.dots.pdf", width = 10, height = 10)

#===========================================================================================================================
# Subset fullloid ligand interactions
full.lig.cells <- levels(exp.m.s.hla.spdb.means$ligand.cells)[grep("mac|Foam", levels(exp.m.s.hla.spdb.means$ligand.cells))]
exp.m.sig.full.lig.hla.spdb.means <- exp.m.s.hla.spdb.means[exp.m.s.hla.spdb.means$ligand.cells %in% full.lig.cells, ]
head(exp.m.sig.full.lig.hla.spdb.means)

# Check dimensions
length(unique(exp.m.sig.full.lig.hla.spdb.means$ligand))
length(unique(exp.m.sig.full.lig.hla.spdb.means$receptor))

# And plot
ggplot(exp.m.sig.full.lig.hla.spdb.means, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("fullloid ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in fullloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("cellphonedb_results/HLA.fullloid.ligands.dots.pdf", width = 10, height = 10)

# Subset a bit more for clarity
summary(exp.m.sig.full.lig.hla.spdb.means$interaction.score)
sub.exp.m.sig.full.lig.hla.spdb.means <- exp.m.sig.full.lig.hla.spdb.means[exp.m.sig.full.lig.hla.spdb.means$interaction.score > summary(exp.m.sig.full.lig.hla.spdb.means$interaction.score)[5],]

summary(summary(exp.m.sig.full.lig.hla.spdb.means$ligand))
sub.exp.m.sig.full.lig.hla.spdb.means <- sub.exp.m.sig.full.lig.hla.spdb.means[sub.exp.m.sig.full.lig.hla.spdb.means$ligand %in% names(summary(sub.exp.m.sig.full.lig.hla.spdb.means$ligand))[summary(summary(sub.exp.m.sig.full.lig.hla.spdb.means$ligand)) < 10],]

# Check dimensions
length(unique(sub.exp.m.sig.full.lig.hla.spdb.means$ligand))
length(unique(sub.exp.m.sig.full.lig.hla.spdb.means$receptor))

# And plot
ggplot(sub.exp.m.sig.full.lig.hla.spdb.means, aes(x = ligand, y= receptor, col = receptor.cells, size = interaction.score, shape = ligand.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("fullloid ligands")) +
  theme_light() +
  ylab("Receptors (in other cells)") +
  xlab("Ligands (in fullloid cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 14),
        aspect.ratio = 1)
ggsave("cellphonedb_results/HLA.fullloid.ligands.subset.dots.pdf", width = 15, height = 15)


#===========================================================================================================================
# Subset fullloid receptor interactions
full.rec.cells <- levels(exp.m.s.hla.spdb.means$receptor.cells)[grep("mac|Foam", levels(exp.m.s.hla.spdb.means$receptor.cells))]
exp.m.sig.full.rec.hla.spdb.means <- exp.m.s.hla.spdb.means[exp.m.s.hla.spdb.means$receptor.cells %in% full.rec.cells, ]
head(exp.m.sig.full.rec.hla.spdb.means)

# Check dimensions
length(unique(exp.m.sig.full.rec.hla.spdb.means$ligand))
length(unique(exp.m.sig.full.rec.hla.spdb.means$receptor))

# And plot
ggplot(exp.m.sig.full.rec.hla.spdb.means, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("fullloid receptors")) +
  theme_light() +
  ylab("Receptors (in fullloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("cellphonedb_results/HLA.fullloid.receptors.dots.pdf", width = 25, height = 25)

# Subset a bit more for clarity
summary(exp.m.sig.full.rec.hla.spdb.means$interaction.score)
sub.exp.m.sig.full.rec.hla.spdb.means <- exp.m.sig.full.rec.hla.spdb.means[exp.m.sig.full.rec.hla.spdb.means$interaction.score > summary(exp.m.sig.full.rec.hla.spdb.means$interaction.score)[5],]

summary(summary(exp.m.sig.full.rec.hla.spdb.means$ligand))
sub.exp.m.sig.full.rec.hla.spdb.means <- sub.exp.m.sig.full.rec.hla.spdb.means[sub.exp.m.sig.full.rec.hla.spdb.means$ligand %in% names(summary(sub.exp.m.sig.full.rec.hla.spdb.means$ligand))[summary(sub.exp.m.sig.full.rec.hla.spdb.means$ligand) < 20],]

# Check dimensions
length(unique(sub.exp.m.sig.full.rec.hla.spdb.means$ligand))
length(unique(sub.exp.m.sig.full.rec.hla.spdb.means$receptor))

# And plot
ggplot(sub.exp.m.sig.full.rec.hla.spdb.means, aes(x = ligand, y= receptor, col = ligand.cells, size = interaction.score, shape = receptor.cells)) +
  geom_point() +
  geom_jitter() +
  scale_fill_manual(values = subclus.colors, aesthetics = "col") +
  ggtitle(paste("fullloid receptors")) +
  theme_light() +
  ylab("Receptors (in fullloid cells)") +
  xlab("Ligands (in other cells)") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        text = element_text(size = 16),
        aspect.ratio = 1)
ggsave("cellphonedb_results/HLA.fullloid.receptors.subset.dots.pdf", width = 15, height = 15)

