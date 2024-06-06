#======================================================================
## Investigate Mosaic Loss of Y
#======================================================================
## Setup
dir.create("mosaicY_monocytes", showWarnings = F)

## Get Y genes
# Extract datatable from Hg db
egCHR    <- toTable(org.Hs.egCHR)
egSYMBOL <- toTable(org.Hs.egSYMBOL)
egCHRSYM <- merge(egCHR, egSYMBOL, by = "gene_id")
head(egCHRSYM)

# Subset The Y genes
egY <- egCHRSYM[egCHRSYM$chromosome == "Y","symbol"]

# Remove lncRNA junk
egY <- egY[grep("LOC|TTTY", egY, invert = T)]

## Make a control set
# Subset all but Y genes
egNOTY <- egCHRSYM[egCHRSYM$chromosome != "Y","symbol"]

# Remove lncRNA junk
egNOTY <- egNOTY[grep("LOC|TTTY|RPL|LINC|DNA|-AS|MIR|RPS|SNOR|orf", egNOTY, invert = T)]
## Keep only genes expressed in the monocytes
# Get the monocytes
mono.pops <- as.vector(unique(Idents(mono.seurat))[grep("Monocytes", unique(Idents(mono.seurat)))])
mono.seurat <- subset(mono.seurat, ident = mono.pops)

# Get the mono expressed gene set
egALL      <- row.names(mono.seurat@assays$RNA@counts)
egALL.expr <- GetAssayData(mono.seurat, slot = "data", assay = "RNA")[egALL,]
egALL      <- egALL[egALL %in% row.names(egALL.expr[apply(egALL.expr, 1, mean) != 0,])]

# Subset the Y genes
egY <- egY[egY %in% egALL]

# Clean up
rm(egALL.expr)

#======================================================================
## Aggregate chr Y gene expression per patient per cell
# Retrieve expression data per cell
egY.expr    <- GetAssayData(mono.seurat, slot = "data", assay = "RNA")[egY,]

# Take the average expression of every cell
egY.expr.mean.cell    <- apply(egY.expr, 2, mean)

## Aggregate it back to patient level
# Get cell names
egY.cells    <- names(egY.expr.mean.cell)

# Get 10X patient names
egY.cells.10X <- egY.cells[grep("PBMC", egY.cells)]
egY.patients  <- paste0("P", substr(egY.cells.10X, nchar(egY.cells.10X), nchar(egY.cells.10X)), sep = "")

## Plot the aggregate gene expression at the patient level
ggboxplot(m.egY.expr.mean.cell, x = "Patient", y = "value", title = "chrY genes average expression per cell per patient", xlab = "Patient", ylab = "Normalized expression", add = c("point", "jitter")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mosaicY_monocytes/chrY_genes_per_patient_per_cell.boxplot.pdf")


#======================================================================
## Determine the ratio of cells with loss of Y
egY.patient.ratio <- list()
for(thePatient in unique(egY.patients)){
  egY.patient.ratio[[thePatient]] <- subset(m.egY.expr.mean.cell, Patient == thePatient)
  p.cells         <- nrow(egY.patient.ratio[[thePatient]])
  p.cells.zero    <- sum(egY.patient.ratio[[thePatient]]$value == 0)
  egY.patient.ratio[[thePatient]] <- p.cells.zero / p.cells
}

# Plot the chr Y loss ratios
m.egY.patient.ratio <- melt(egY.patient.ratio)
colnames(m.egY.patient.ratio) <- c("loss.ratio", "Patient")
m.egY.patient.ratio <- m.egY.patient.ratio[order(m.egY.patient.ratio$loss.ratio, decreasing = T),]
ggbarplot(m.egY.patient.ratio, x = "Patient", y = "loss.ratio", title = "chrY genes ratio of cells with 0 reads on Y per (male) patient", xlab = "Patient", ylab = "Loss Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mosaicY_monocytes/chrY_genes_ratio_per_patient.barplot.pdf")

#======================================================================
## Now explore the ratio's between cells with presumed loss of Y versus not as opposed to the controls
# Determine the ratio of cells with loss of Y
egY.patient.ratio <- list()
egY.patient.cells <- list()
egY.zero.cells    <- list()
egY.notzero.cells <- list()
for(thePatient in unique(egY.patients)){
  egY.patient.cells[[thePatient]] <- subset(m.egY.expr.mean.cell, Patient == thePatient)
  p.cells         <- nrow(egY.patient.cells[[thePatient]])
  p.cells.zero    <- sum(egY.patient.cells[[thePatient]]$value == 0)
  egY.zero.cells[[thePatient]]    <- row.names(egY.patient.cells[[thePatient]][egY.patient.cells[[thePatient]]$value == 0,])
  egY.notzero.cells[[thePatient]] <- row.names(egY.patient.cells[[thePatient]][egY.patient.cells[[thePatient]]$value != 0,])
  egY.patient.ratio[[thePatient]] <- p.cells.zero / p.cells
}

# Plot the average expression in zero Y cells
avg.exp.all_but_Y <- as.data.frame(GetAssayData(mono.seurat, slot = "data", assay = "RNA")[,unlist(egY.zero.cells)])
avg.exp.all_but_Y <- avg.exp.all_but_Y[which(!row.names(avg.exp.all_but_Y) %in% egY),]
pdf("mosaicY_monocytes/not_on_chrY_Genes_avgexp_in_no_chrY_cells.hist.pdf")
hist(apply(avg.exp.all_but_Y, 2, mean), main = "Average expression of all but chrY genes in cells with no chr Y expression", xlab = "Normalized expression")
dev.off()

# Plot the number of detected transcripts in zero Y cells
VlnPlot(object = subset(mono.seurat, cells = unlist(egY.zero.cells)), group.by = "Patient", features = c("nFeature_RNA","nCount_RNA"), y.max = 15000) + 
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold"), 
        aspect.ratio     = 1
  )
ggsave("mosaicY_monocytes/featureCount.chrY_Zero cells.pdf")

# Plot the average expression in nonzero Y cells
avg.exp.nonzeroY <- as.data.frame(GetAssayData(mono.seurat, slot = "data", assay = "RNA")[,unlist(egY.notzero.cells)])
avg.exp.nonzeroY <- avg.exp.nonzeroY[which(!row.names(avg.exp.nonzeroY) %in% egY),]
pdf("mosaicY_monocytes/not_on_chrY_Genes_avgexp_in_notzero_chrY_cells.hist.pdf")
hist(apply(avg.exp.nonzeroY, 2, mean), main = "Average expression of all but chrY genes in cells with no chr Y expression", xlab = "Normalized expression")
dev.off()

# Plot the number of detected transcripts in nonzero Y cells
VlnPlot(object = subset(mono.seurat, cells = unlist(egY.notzero.cells)), group.by = "Patient", features = c("nFeature_RNA","nCount_RNA"), y.max = 15000) + 
  theme_pubr(base_size = 14, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold"), 
        aspect.ratio     = 1
  )
ggsave("mosaicY_monocytes/featureCount.chrY_notZero cells.pdf")

## Plot the absolute number of zero and non zero Y cells per patient
# Prep the data
egY.cells.abs         <- data.frame(Zero = sort(unlist(lapply(egY.zero.cells,length)), decreasing = T), Not.zero = sort(unlist(lapply(egY.notzero.cells,length)), decreasing = T))
egY.cells.abs$Patient <- row.names(egY.cells.abs)
m.egY.cells.abs       <- melt(egY.cells.abs)
colnames(m.egY.cells.abs) <- c("Patient", "chrY.genes", "cell.count")

# Plot
ggbarplot(m.egY.cells.abs[m.egY.cells.abs$Patient %in% c("P1", "P2", "P3"),], x = "Patient", y = "cell.count", fill = "chrY.genes", title = "Cells with or without chrY gene expression", xlab = "Patient", ylab = "Cell numbers", palette = "Paired", position = position_dodge2()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mosaicY_monocytes/chrY_genes_ratio_per_patient.10X.barplot.pdf")


#======================================================================
## Bootstrap the analysis of the CTRL genes
CTRL.boots.expr.mean.cell <- list()
CTRL.boots.genes          <- list()
for(n in 1:1000){
  cat("Bootstrapping: ", n, "\n")
  
  # Subset all but Y genes
  egNOTY <- egCHRSYM[egCHRSYM$chromosome != "Y","symbol"]
  
  # Remove lncRNA junk
  egNOTY <- egNOTY[grep("LOC|TTTY|RPL|LINC|DNA|-AS|MIR|RPS|SNOR|orf", egNOTY, invert = T)]
  
  # Subset the control genes
  egNOTY <- egNOTY[egNOTY %in% egALL]
  
  # Keep a random fraction of control genes to match the number of Y genes
  CTRL.boots.genes[[n]] <- sample(egNOTY, size = length(egY), replace = F)
  
  #======================================================================
  ## Aggregate gene expression per patient per cell
  # Retrieve expression data per cell
  tmp.egNOTY.expr <- as.data.frame(GetAssayData(mono.seurat, slot = "data", assay = "RNA")[CTRL.boots.genes[[n]],])
  
  # Take the average expression of every cell
  tmp.egNOTY.expr <- apply(tmp.egNOTY.expr, 2, mean)
  
  ## Aggregate it back to patient level
  # Get cell names
  egNOTY.cells <- names(tmp.egNOTY.expr)
  
  # Get 10X patient names
  egNOTY.patients <- egY.patients
 
  # Melt into a frame with patient info
  tmp.m.expr         <- melt(tmp.egNOTY.expr)
  tmp.m.expr$Patient <- egNOTY.patients
  
  # Delete the females
  CTRL.boots.expr.mean.cell[[n]] <- tmp.m.expr[tmp.m.expr$Patient %in% egY.patients,]
}

#======================================================================
## Now explore the ratio's between cells with presumed loss of Y versus not as opposed to the controls
CTRL.boots.patient.ratio <- list()
m.patient.ratio          <- list()
Y.vs.CTRL.pval           <- list()
ordered.patient.ratio    <- list()
for(n in 1:1000){
  cat("Bootstrapping: ", n, "\n")
  
  # Determine the ratio of cells with no reads on CTRL genes
  tmp.CTRL <- list()
  for(thePatient in unique(egY.patients)){
    tmp.CTRL[[thePatient]] <- subset(CTRL.boots.expr.mean.cell[[n]], Patient == thePatient)
    p.cells         <- nrow(tmp.CTRL[[thePatient]])
    p.cells.zero    <- sum(tmp.CTRL[[thePatient]]$value == 0)
    tmp.CTRL[[thePatient]] <- p.cells.zero / p.cells
  }
  
  # Tidy up
  CTRL.boots.patient.ratio[[n]] <- melt(tmp.CTRL)
  colnames(CTRL.boots.patient.ratio[[n]]) <- c("loss.ratio", "Patient")
  
  # Merge the two sets for easier comparison
  m.patient.ratio[[n]] <- m.egY.patient.ratio
  m.patient.ratio[[n]]$CTRL <- CTRL.boots.patient.ratio[[n]]$loss.ratio
  colnames(m.patient.ratio[[n]]) <- c("chrY", "Patient", "CTRL")
  m.patient.ratio[[n]] <- melt(m.patient.ratio[[n]])
  colnames(m.patient.ratio[[n]]) <- c("Patient", "gene.set", "value")
  
  # Check for a statistical difference
  Y.vs.CTRL.pval[[n]] <- wilcox.test(m.egY.patient.ratio$loss.ratio, CTRL.boots.patient.ratio[[n]]$loss.ratio, paired = T, alternative = "g")$p.value
  
  # Order on  ratio diff
  ordered.patient.ratio[[n]]            <- m.egY.patient.ratio
  ordered.patient.ratio[[n]]$diff       <- m.egY.patient.ratio$loss.ratio - CTRL.boots.patient.ratio[[n]]$loss.ratio
  ordered.patient.ratio[[n]]$loss.ratio <- NULL
  ordered.patient.ratio[[n]]            <- ordered.patient.ratio[[n]][order(ordered.patient.ratio[[n]]$diff, decreasing = F),]
}

median(unlist(Y.vs.CTRL.pval))
Y.vs.CTRL.padj <- p.adjust(unlist(Y.vs.CTRL.pval), method = "BH")
sum(Y.vs.CTRL.padj < 0.05)

#======================================================================
## Take the median and the SD of the CTRL bootstraps
# Populate a dataframe with all the bootstrap values per patient
CTRL.boots.patient.ratio.df <- m.patient.ratio[[1]][m.patient.ratio[[n]]$gene.set == "CTRL", c("Patient", "value")]
colnames(CTRL.boots.patient.ratio.df) <- c("Patient", "bs1")
for(n in 2:1000){
  cat("Bootstrap", n, "\n")
  cn <- colnames(CTRL.boots.patient.ratio.df)
  CTRL.boots.patient.ratio.df <- merge(CTRL.boots.patient.ratio.df, m.patient.ratio[[n]][m.patient.ratio[[n]]$gene.set == "CTRL", c("Patient", "value")], by = "Patient")
  colnames(CTRL.boots.patient.ratio.df) <- c(cn, paste("bs", n, sep = ""))
}

# Tidy up the df
rownames(CTRL.boots.patient.ratio.df) <- CTRL.boots.patient.ratio.df$Patient
CTRL.boots.patient.ratio.df$Patient   <- NULL

# Take the median value per patient
CTRL.boots.patient.ratio.median <- data.frame(Patient  = names(apply(CTRL.boots.patient.ratio.df,1, median)), 
                                              gene.set = rep("CTRL", length(apply(CTRL.boots.patient.ratio.df,1, median))),
                                              value    = apply(CTRL.boots.patient.ratio.df,1, median))
row.names(CTRL.boots.patient.ratio.median) <- NULL

# Add back to a molten df with the chrY observations per patient
tmp.patient.ratio <- m.patient.ratio[[1]][m.patient.ratio[[n]]$gene.set == "chrY", ]
tmp.patient.ratio <- tmp.patient.ratio[order(tmp.patient.ratio$value, decreasing = T),]
m.patient.ratio.bs.median <- rbind(tmp.patient.ratio, CTRL.boots.patient.ratio.median)

# Calculate significance
Y.vs.CTRL.pval.bs.median <- wilcox.test(m.egY.patient.ratio$loss.ratio, CTRL.boots.patient.ratio.median$value, paired = T, alternative = "g")$p.value
median(unlist(Y.vs.CTRL.padj))

# Take the SD
CTRL.boots.patient.ratio.sd <- data.frame(row.names = row.names(CTRL.boots.patient.ratio.df), SD = rep(0, nrow(CTRL.boots.patient.ratio.df)))
for(thePatient in row.names(CTRL.boots.patient.ratio.df)){
  CTRL.boots.patient.ratio.sd[thePatient,"SD"] <- sd(CTRL.boots.patient.ratio.df[thePatient,])
}

# Plot both chrY and CTRL sets together
ggbarplot(m.patient.ratio.bs.median, x = "Patient", y = "value", title = "chrY genes ratio of cells with 0 reads on Y per (male) patient", xlab = "Patient", ylab = "Zero Cell Ratio", fill = "gene.set", palette = "Paired", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = mean(m.patient.ratio.bs.median[m.patient.ratio.bs.median$gene.set == "CTRL", "value"]), color = brewer.pal(3, "Paired")[2]) +
  geom_hline(yintercept = mean(m.patient.ratio.bs.median[m.patient.ratio.bs.median$gene.set == "chrY", "value"]), color = brewer.pal(3, "Paired")[1]) +
  annotate("text", x = -Inf, y = Inf, label = paste("One-sided two-sample paired rank sum test p-value: ", round(Y.vs.CTRL.pval.bs.median, digits = 5), sep = ""), hjust = -0.1, vjust = 1)

ggsave("mosaicY_monocytes/genes_ratio_per_patient_bootstrapped.barplot.pdf")

#======================================================================
## Check robustness per patient
# Initialize dfs
robust.boots.under.df <- data.frame(row.names = ordered.patient.ratio[[n]][, "Patient"])
robust.boots.equal.df <- data.frame(row.names = ordered.patient.ratio[[n]][, "Patient"])
robust.boots.over.df  <- data.frame(row.names = ordered.patient.ratio[[n]][, "Patient"])

# Count occurence of directionality
for(n in 1:1000){
  cat("Bootstrap", n, "\n")
  robust.boots.under.df[ordered.patient.ratio[[n]][ordered.patient.ratio[[n]]$diff <  0, "Patient"], paste("bs", n, sep = "")] <- 1
  robust.boots.equal.df[ordered.patient.ratio[[n]][ordered.patient.ratio[[n]]$diff == 0, "Patient"], paste("bs", n, sep = "")] <- 1
  robust.boots.over.df[ ordered.patient.ratio[[n]][ordered.patient.ratio[[n]]$diff >  0, "Patient"], paste("bs", n, sep = "")] <- 1
}

# Tidy up
robust.boots.under.df[is.na(robust.boots.under.df)] <- 0
robust.boots.equal.df[is.na(robust.boots.equal.df)] <- 0
robust.boots.over.df[ is.na(robust.boots.over.df)]  <- 0

# Sum the occurence
robust.boots.under <- apply(robust.boots.under.df, 1, sum)
robust.boots.equal <- apply(robust.boots.equal.df, 1, sum)
robust.boots.over  <- apply(robust.boots.over.df,  1, sum)

# Asymptotic Model fit
set.seed(666)
X <- 1:length(robust.boots.over)
Y <- sort(robust.boots.over) / 10
model <- drm(Y ~ X, fct = DRC.asymReg())
pdf("mosaicY_monocytes/bootstrap_fit_cutoff.scatter.pdf")
plot(model, log="", main = "Inflection point of fitted curve on robustness of bootstraps", xlab = "Patient", ylab = "% of bootstraps", )
find_curve_elbow(data_frame = data.frame(X, Y), export_type = "all")
curve.fit.cutoff <- find_curve_elbow(data_frame = data.frame(X, Y), export_type = "all")["Y"]
abline(v = find_curve_elbow(data_frame = data.frame(X, Y)), col = "red")
dev.off()

# Populate a df
robust.boots.df <- data.frame(under = robust.boots.under, equal = robust.boots.equal, over = robust.boots.over)
robust.boots.df <- robust.boots.df[order(robust.boots.df$over),]

# Melt it for plotting
m.robust.boots.df         <- robust.boots.df
m.robust.boots.df$Patient <- row.names(m.robust.boots.df)
m.robust.boots.df         <- melt(m.robust.boots.df)
m.robust.boots.df$value   <- m.robust.boots.df$value / 10
colnames(m.robust.boots.df) <- c("Patient", "Direction", "Percentage")

ggbarplot(m.robust.boots.df, x = "Patient", y = "Percentage", fill = "Direction", title = "Robustness of bootstrap representation", xlab = "Patient", ylab = "% of bootstraps", palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = curve.fit.cutoff, color = "gold", linewidth = 0.75)
ggsave("mosaicY_monocytes/bootstrap_directionality.barplot.pdf")

# Construct ratio diff df
ordered.patient.ratio.bs.median            <- m.egY.patient.ratio
ordered.patient.ratio.bs.median$diff       <- m.egY.patient.ratio$loss.ratio - CTRL.boots.patient.ratio.median$value
ordered.patient.ratio.bs.median$loss.ratio <- NULL
ordered.patient.ratio.bs.median            <- merge(ordered.patient.ratio.bs.median, CTRL.boots.patient.ratio.sd, by.y = 0, by.x = "Patient")
ordered.patient.ratio.bs.median            <- ordered.patient.ratio.bs.median[order(ordered.patient.ratio.bs.median$diff, decreasing = F),]

# Add in robustness of overrepresentation
m.over_cutoff.robust.boots <- m.robust.boots.df[m.robust.boots.df$Direction == "over",]
m.over_cutoff.robust.boots <- m.over_cutoff.robust.boots[m.over_cutoff.robust.boots[, "Percentage"] >= curve.fit.cutoff, "Patient"]
ordered.patient.ratio.bs.median$Robust <- "No"
ordered.patient.ratio.bs.median[ordered.patient.ratio.bs.median$Patient %in% m.over_cutoff.robust.boots, "Robust"] <- "Yes"

# Plot ratio diff with robustness cutoff
ggbarplot(ordered.patient.ratio.bs.median, x = "Patient", y = "diff", fill = "Robust", title = "Ratio differential of chrY vs. CTRL genes", xlab = "Patient", ylab = 'Difference in "zero cell" ratio (chrY - CTRL)', palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_errorbar(aes(ymin=diff-SD, ymax=diff+SD), width=.2,
                position=position_dodge(.9))
ggsave("mosaicY_monocytes/genes_ratio_diff_per_patient_bootstrapped.barplot.pdf")


#======================================================================
## Feed mosaic Y metadata back to the Seurat Object
## Patient level






