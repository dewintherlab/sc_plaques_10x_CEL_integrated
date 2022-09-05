#======================================================================
## Make pretty plots of the initital populations
#======================================================================
## 10X data 
## First need to rename the populations
celltypes.10X <- refined.patients.seurat$predicted.celltype.l2
unique(celltypes.10X)

# T cells
celltypes.10X[grep("CD4", celltypes.10X)]  <- "T Cells"
celltypes.10X[grep("CD8", celltypes.10X)]  <- "T Cells"
celltypes.10X[grep("Treg", celltypes.10X)] <- "T Cells"
celltypes.10X[grep("dnT", celltypes.10X)]  <- "T Cells"
celltypes.10X[grep("gdT", celltypes.10X)]  <- "T Cells"
celltypes.10X[grep("MAIT", celltypes.10X)] <- "T Cells"

# B cells
celltypes.10X[grep("B ", celltypes.10X)] <- "B Cells"

# NK Cells
celltypes.10X[grep("NK", celltypes.10X)] <- "NK Cells"

# Myeloid cells
celltypes.10X[grep("CD14", celltypes.10X)] <- "Myeloid Cells"
celltypes.10X[grep("CD16", celltypes.10X)] <- "Myeloid Cells"
celltypes.10X[grep("DC", celltypes.10X)]   <- "Myeloid Cells"

# Junk
celltypes.10X[grep("Eryth", celltypes.10X)]    <- "Erythroid Cells"
celltypes.10X[grep("ILC", celltypes.10X)]      <- "Innate Lymphoid Cells"
celltypes.10X[grep("Platelet", celltypes.10X)] <- "Platelets"

# Add back to the seurat object
refined.patients.seurat$celltypes.10X <- celltypes.10X

# Define some colors for these pops
unique(celltypes.10X)

# T
celltypes.10X.cols <- "orange3"

# B
celltypes.10X.cols <- c(celltypes.10X.cols, "purple")

# Plasmablast
celltypes.10X.cols <- c(celltypes.10X.cols, "azure2")

# NK
celltypes.10X.cols <- c(celltypes.10X.cols, "mediumslateblue")

# Myeloid
celltypes.10X.cols <- c(celltypes.10X.cols, "dodgerblue2")

# HSPC
celltypes.10X.cols <- c(celltypes.10X.cols, "azure3")

# Ery
celltypes.10X.cols <- c(celltypes.10X.cols, "red4")

# Doublet
celltypes.10X.cols <- c(celltypes.10X.cols, "azure")

# Platelets
celltypes.10X.cols <- c(celltypes.10X.cols, "azure4")

# ILC
celltypes.10X.cols <- c(celltypes.10X.cols, "bisque1")

# Name the pops
names(celltypes.10X.cols) <- unique(celltypes.10X)


## And draw the plot
customUMAP(object = refined.patients.seurat, 
           pt.size = 0.5,
           group.by = "celltypes.10X",
           title = "10X patients", cols = celltypes.10X.cols,
           shuffle = T, 
           legend.pos = "right", 
           seed = 666, plot.width = 10, plot.height = 10,
           file.name = "integration_QC/Integration UMAP 10X.pdf")


##====================================================================================
## CEL-seq data
## First need to rename the populations
celltypes.43p <- as.vector(Idents(full.43p.seurat))
unique(celltypes.43p)

# T cells
celltypes.43p[grep("T ", celltypes.43p)] <- "T Cells"
celltypes.43p[grep("FOXP3", celltypes.43p)]  <- "T Cells"

# B cells
celltypes.43p[grep("CD79", celltypes.43p)] <- "B Cells"

# NK Cells
celltypes.43p[grep("NK", celltypes.43p)] <- "NK Cells"

# Myeloid cells
celltypes.43p[grep("CD68", celltypes.43p)] <- "Myeloid Cells"

# SMC
celltypes.43p[grep("ACTA", celltypes.43p)] <- "Smooth Muscle Cells"

# Endothelial cells
celltypes.43p[grep("Endo", celltypes.43p)] <- "Endothelial Cells"


# Add back to the seurat object
full.43p.seurat$celltypes.43p <- celltypes.43p

# Define some colors for these pops
unique(celltypes.43p)

# T
celltypes.43p.cols <- "orange3"

# Endothelial
celltypes.43p.cols <- c(celltypes.43p.cols, "chartreuse4")

# NK
celltypes.43p.cols <- c(celltypes.43p.cols, "mediumslateblue")

# Myeloid
celltypes.43p.cols <- c(celltypes.43p.cols, "dodgerblue2")

# SMC
celltypes.43p.cols <- c(celltypes.43p.cols, "chartreuse2")

# B
celltypes.43p.cols <- c(celltypes.43p.cols, "purple")


# Name the pops
names(celltypes.43p.cols) <- unique(celltypes.43p)


## And draw the plot
customUMAP(object = full.43p.seurat, 
           pt.size = 0.5,
           group.by = "celltypes.43p",
           title = "43p patients", cols = celltypes.43p.cols,
           shuffle = T, 
           legend.pos = "right", 
           seed = 666, plot.width = 10, plot.height = 10,
           file.name = "integration_QC/Integration UMAP CEL-seq.pdf")



