#======================================================================
## Investigate Mosaic Loss of Y
#======================================================================
## Setup
dir.create("mosaicY_macrophages", showWarnings = F)

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

## Keep only genes expressed in the myeloid cells
# Get the myeloid expressed gene set
egALL      <- row.names(final.pop.call.from_full.integrated.mac.seurat@assays$RNA@counts)
egALL.expr <- GetAssayData(final.pop.call.from_full.integrated.mac.seurat, slot = "data", assay = "RNA")[egALL,]
egALL      <- egALL[egALL %in% row.names(egALL.expr[apply(egALL.expr, 1, mean) != 0,])]

# Subset the Y genes
egY <- egY[egY %in% egALL]

# Clean up
rm(egALL.expr)

#======================================================================
## Aggregate chr Y gene expression per patient per cell
# Retrieve expression data per cell
egY.expr    <- GetAssayData(final.pop.call.from_full.integrated.mac.seurat, slot = "data", assay = "RNA")[egY,]

# Take the average expression of every cell
egY.expr.mean.cell    <- apply(egY.expr, 2, mean)

## Aggregate it back to patient level
# Get cell names
egY.cells    <- names(egY.expr.mean.cell)

# Get 10X patient names
egY.cells.10X    <- egY.cells[grep("plaque|PBMC", egY.cells)]
egY.patients.10X <- paste0("P", substr(egY.cells.10X, nchar(egY.cells.10X), nchar(egY.cells.10X)), sep = "")

# Get CEL-seq patient names
egY.cells.CS    <- egY.cells[grep("plaque|PBMC", egY.cells, invert = T)]
egY.patients.CS <- paste0("P", substr(egY.cells.CS, 1, 4), sep = "")

# Merge all patient names
egY.patients    <- c(egY.patients.10X, egY.patients.CS)

# Melt into a frame with patient info
m.egY.expr.mean.cell <- melt(egY.expr.mean.cell)
m.egY.expr.mean.cell$Patient <- egY.patients


## Remove female patients
# Get male patients
male.cells <- row.names(final.pop.call.from_full.integrated.mac.seurat@meta.data[final.pop.call.from_full.integrated.mac.seurat@meta.data$Sex == "male",])

# Get 10X patient names
male.cells.10X    <- male.cells[grep("plaque|PBMC", male.cells)]
male.patients.10X <- paste0("P", substr(male.cells.10X, nchar(male.cells.10X), nchar(male.cells.10X)), sep = "")

# Get CEL-seq patient names
male.cells.CS    <- male.cells[grep("plaque|PBMC", male.cells, invert = T)]
male.patients.CS <- paste0("P", substr(male.cells.CS, 1, 4), sep = "")

# Merge all patient names
male.patients    <- unique(c(male.patients.10X, male.patients.CS))
male.patients

# Delete the females
m.egY.expr.mean.cell    <- m.egY.expr.mean.cell[m.egY.expr.mean.cell$Patient %in% male.patients,]

## Plot the aggregate gene expression at the patient level
ggboxplot(m.egY.expr.mean.cell, x = "Patient", y = "value", title = "chrY genes average expression per cell per patient", xlab = "Patient", ylab = "Normalized expression", add = c("point", "jitter")) +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mosaicY_macrophages/chrY_genes_per_patient_per_cell.boxplot.pdf")


#======================================================================
## Determine the ratio of cells with loss of Y
egY.patient.ratio <- list()
for(thePatient in unique(male.patients)){
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
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mosaicY_macrophages/chrY_genes_ratio_per_patient.barplot.pdf")

#======================================================================
## Now explore the ratio's between cells with presumed loss of Y versus not as opposed to the controls
# Determine the ratio of cells with loss of Y
egY.patient.ratio <- list()
egY.patient.cells <- list()
egY.zero.cells    <- list()
egY.notzero.cells <- list()
for(thePatient in unique(male.patients)){
  egY.patient.cells[[thePatient]] <- subset(m.egY.expr.mean.cell, Patient == thePatient)
  p.cells         <- nrow(egY.patient.cells[[thePatient]])
  p.cells.zero    <- sum(egY.patient.cells[[thePatient]]$value == 0)
  egY.zero.cells[[thePatient]]    <- row.names(egY.patient.cells[[thePatient]][egY.patient.cells[[thePatient]]$value == 0,])
  egY.notzero.cells[[thePatient]] <- row.names(egY.patient.cells[[thePatient]][egY.patient.cells[[thePatient]]$value != 0,])
  egY.patient.ratio[[thePatient]] <- p.cells.zero / p.cells
}

# Plot the average expression in zero Y cells
avg.exp.all_but_Y <- as.data.frame(GetAssayData(final.pop.call.from_full.integrated.mac.seurat, slot = "data", assay = "RNA")[,unlist(egY.zero.cells)])
avg.exp.all_but_Y <- avg.exp.all_but_Y[which(!row.names(avg.exp.all_but_Y) %in% egY),]
pdf("mosaicY_macrophages/not_on_chrY_Genes_avgexp_in_no_chrY_cells.hist.pdf")
hist(apply(avg.exp.all_but_Y, 2, mean), main = "Average expression of all but chrY genes in cells with no chr Y expression", xlab = "Normalized expression")
dev.off()

# Plot the number of detected transcripts in zero Y cells
VlnPlot(object = subset(final.pop.call.from_full.integrated.mac.seurat, cells = unlist(egY.zero.cells)), group.by = "Patient", features = c("nFeature_RNA","nCount_RNA"), y.max = 15000) + 
  theme_pubr(base_size = 16, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold"), 
        aspect.ratio     = 1
  )
ggsave("mosaicY_macrophages/featureCount.chrY_Zero cells.pdf")

# Plot the average expression in nonzero Y cells
avg.exp.nonzeroY <- as.data.frame(GetAssayData(final.pop.call.from_full.integrated.mac.seurat, slot = "data", assay = "RNA")[,unlist(egY.notzero.cells)])
avg.exp.nonzeroY <- avg.exp.nonzeroY[which(!row.names(avg.exp.nonzeroY) %in% egY),]
pdf("mosaicY_macrophages/not_on_chrY_Genes_avgexp_in_notzero_chrY_cells.hist.pdf")
hist(apply(avg.exp.nonzeroY, 2, mean), main = "Average expression of all but chrY genes in cells with no chr Y expression", xlab = "Normalized expression")
dev.off()

# Plot the number of detected transcripts in nonzero Y cells
VlnPlot(object = subset(final.pop.call.from_full.integrated.mac.seurat, cells = unlist(egY.notzero.cells)), group.by = "Patient", features = c("nFeature_RNA","nCount_RNA"), y.max = 15000) + 
  theme_pubr(base_size = 16, x.text.angle = 45, legend = "right") +
  theme(panel.background = element_blank(),
        title            = element_text(size = 16, face = "bold"), 
        aspect.ratio     = 1
  )
ggsave("mosaicY_macrophages/featureCount.chrY_notZero cells.pdf")

## Plot the absolute number of zero and non zero Y cells per patient
# Prep the data
egY.cells.abs         <- data.frame(Zero = sort(unlist(lapply(egY.zero.cells,length)), decreasing = T), Not.zero = sort(unlist(lapply(egY.notzero.cells,length)), decreasing = T))
egY.cells.abs$Patient <- row.names(egY.cells.abs)
m.egY.cells.abs       <- melt(egY.cells.abs)
colnames(m.egY.cells.abs) <- c("Patient", "chrY.genes", "cell.count")

# Plot, but split off the 10X samples for visibility
ggbarplot(m.egY.cells.abs[m.egY.cells.abs$Patient %in% c("P1", "P2", "P3"),], x = "Patient", y = "cell.count", fill = "chrY.genes", title = "Cells with or without chrY gene expression", xlab = "Patient", ylab = "Cell numbers", palette = "Paired", position = position_dodge2()) +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mosaicY_macrophages/chrY_genes_ratio_per_patient.10X.barplot.pdf")

ggbarplot(m.egY.cells.abs[!m.egY.cells.abs$Patient %in% c("P1", "P2", "P3"),], x = "Patient", y = "cell.count", fill = "chrY.genes", title = "Cells with or without chrY gene expression", xlab = "Patient", ylab = "Cell numbers", palette = "Paired", position = position_dodge2()) +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mosaicY_macrophages/chrY_genes_ratio_per_patient.CS.barplot.pdf")

sum(m.egY.cells.abs[!m.egY.cells.abs$Patient %in% c("P1", "P2", "P3"),"cell.count"])

#======================================================================
# Final check: control genes from the same cells as the zero Y genes
# Bootstrap the control gens
same_cell.CTRL.boots.expr.mean.cell <- list()
same_cell.CTRL.boots.genes          <- list()
for(n in 1:1000){
  cat("Bootstrapping: ", n, "\n")
  
  # Subset all but Y genes
  egNOTY <- egCHRSYM[egCHRSYM$chromosome != "Y","symbol"]
  
  # Remove lncRNA junk
  egNOTY <- egNOTY[grep("LOC|TTTY|RPL|LINC|DNA|-AS|MIR|RPS|SNOR|orf", egNOTY, invert = T)]
  
  # Subset the control genes
  egNOTY <- egNOTY[egNOTY %in% egALL]
  
  # Keep a random fraction of control genes to match the number of Y genes
  same_cell.CTRL.boots.genes[[n]] <- sample(egNOTY, size = length(egY), replace = F)
  
  #======================================================================
  ## Aggregate gene expression per patient per cell
  # Retrieve expression data per cell
  tmp.egNOTY.expr <- as.data.frame(GetAssayData(final.pop.call.from_full.integrated.mac.seurat, slot = "data", assay = "RNA")[same_cell.CTRL.boots.genes[[n]],unlist(egY.zero.cells)])
  
  # Take the average expression of every cell
  tmp.egNOTY.expr <- apply(tmp.egNOTY.expr, 2, mean)
  
  ## Aggregate it back to patient level
  # Get cell names
  egNOTY.cells <- names(tmp.egNOTY.expr)
  
  # Get 10X patient names
  egNOTY.cells.10X    <- egNOTY.cells[grep("plaque", egNOTY.cells)]
  egNOTY.patients.10X <- paste0("P", substr(egNOTY.cells.10X, nchar(egNOTY.cells.10X), nchar(egNOTY.cells.10X)), sep = "")
  
  # Get CEL-seq patient names
  egNOTY.cells.CS    <- egNOTY.cells[grep("plaque", egNOTY.cells, invert = T)]
  egNOTY.patients.CS <- paste0("P", substr(egNOTY.cells.CS, 1, 4), sep = "")
  
  # Merge all patient names
  egNOTY.patients <- c(egNOTY.patients.10X, egNOTY.patients.CS)
  
  # Melt into a frame with patient info
  tmp.m.expr         <- melt(tmp.egNOTY.expr)
  tmp.m.expr$Patient <- egNOTY.patients
  
  # Delete the females
  same_cell.CTRL.boots.expr.mean.cell[[n]] <- tmp.m.expr[tmp.m.expr$Patient %in% male.patients,]
}

#======================================================================
## Now explore the ratio's between cells with presumed loss of Y versus not as opposed to the controls
same_cell.CTRL.boots.patient.ratio <- list()
m.patient.ratio          <- list()
Y.vs.CTRL.pval           <- list()
ordered.patient.ratio    <- list()
for(n in 1:1000){
  cat("Bootstrapping: ", n, "\n")
  
  # Determine the ratio of cells with no reads on CTRL genes
  tmp.CTRL <- list()
  for(thePatient in unique(male.patients)){
    tmp.CTRL[[thePatient]] <- subset(same_cell.CTRL.boots.expr.mean.cell[[n]], Patient == thePatient)
    if(nrow(tmp.CTRL[[thePatient]]) != 0){
      p.cells         <- nrow(tmp.CTRL[[thePatient]])
      p.cells.zero    <- sum(tmp.CTRL[[thePatient]]$value == 0)
      tmp.CTRL[[thePatient]] <- p.cells.zero / p.cells
    } else{
      tmp.CTRL[[thePatient]] <- 0
    }
  }
  
  # Tidy up
  same_cell.CTRL.boots.patient.ratio[[n]] <- melt(tmp.CTRL)
  colnames(same_cell.CTRL.boots.patient.ratio[[n]]) <- c("loss.ratio", "Patient")
  
  # Merge the two sets for easier comparison
  m.patient.ratio[[n]] <- m.egY.patient.ratio
  m.patient.ratio[[n]]$CTRL <- same_cell.CTRL.boots.patient.ratio[[n]]$loss.ratio
  colnames(m.patient.ratio[[n]]) <- c("chrY", "Patient", "CTRL")
  m.patient.ratio[[n]] <- melt(m.patient.ratio[[n]])
  colnames(m.patient.ratio[[n]]) <- c("Patient", "gene.set", "value")
  
  # Check for a statistical difference
  Y.vs.CTRL.pval[[n]] <- wilcox.test(m.egY.patient.ratio$loss.ratio, same_cell.CTRL.boots.patient.ratio[[n]]$loss.ratio, paired = T, alternative = "g")$p.value
  
  # Order on  ratio diff
  ordered.patient.ratio[[n]]            <- m.egY.patient.ratio
  ordered.patient.ratio[[n]]$diff       <- m.egY.patient.ratio$loss.ratio - same_cell.CTRL.boots.patient.ratio[[n]]$loss.ratio
  ordered.patient.ratio[[n]]$loss.ratio <- NULL
  ordered.patient.ratio[[n]]            <- ordered.patient.ratio[[n]][order(ordered.patient.ratio[[n]]$diff, decreasing = F),]
}

median(unlist(Y.vs.CTRL.pval))
Y.vs.CTRL.padj <- p.adjust(unlist(Y.vs.CTRL.pval), method = "BH")
sum(Y.vs.CTRL.padj < 0.05)

#======================================================================
## Take the median and the SD of the CTRL bootstraps
# Populate a dataframe with all the bootstrap values per patient
same_cell.CTRL.boots.patient.ratio.df <- m.patient.ratio[[1]][m.patient.ratio[[n]]$gene.set == "CTRL", c("Patient", "value")]
colnames(same_cell.CTRL.boots.patient.ratio.df) <- c("Patient", "bs1")
for(n in 2:1000){
  cat("Bootstrap", n, "\n")
  cn <- colnames(same_cell.CTRL.boots.patient.ratio.df)
  same_cell.CTRL.boots.patient.ratio.df <- merge(same_cell.CTRL.boots.patient.ratio.df, m.patient.ratio[[n]][m.patient.ratio[[n]]$gene.set == "CTRL", c("Patient", "value")], by = "Patient")
  colnames(same_cell.CTRL.boots.patient.ratio.df) <- c(cn, paste("bs", n, sep = ""))
}

# Tidy up the df
rownames(same_cell.CTRL.boots.patient.ratio.df) <- same_cell.CTRL.boots.patient.ratio.df$Patient
same_cell.CTRL.boots.patient.ratio.df$Patient   <- NULL

# Take the median value per patient
same_cell.CTRL.boots.patient.ratio.median <- data.frame(Patient  = names(apply(same_cell.CTRL.boots.patient.ratio.df,1, median)),
                                                        gene.set = rep("CTRL", length(apply(same_cell.CTRL.boots.patient.ratio.df,1, median))),
                                                        value    = apply(same_cell.CTRL.boots.patient.ratio.df,1, median))
row.names(same_cell.CTRL.boots.patient.ratio.median) <- NULL

# Add back to a molten df with the chrY observations per patient
tmp.patient.ratio <- m.patient.ratio[[1]][m.patient.ratio[[n]]$gene.set == "chrY", ]
tmp.patient.ratio <- tmp.patient.ratio[order(tmp.patient.ratio$value, decreasing = T),]
m.patient.ratio.bs.median <- rbind(tmp.patient.ratio, same_cell.CTRL.boots.patient.ratio.median)

# Calculate significance
Y.vs.CTRL.pval.bs.median <- wilcox.test(m.egY.patient.ratio$loss.ratio, same_cell.CTRL.boots.patient.ratio.median$value, paired = T, alternative = "g")$p.value
median(unlist(Y.vs.CTRL.padj))

# Take the SD
same_cell.CTRL.boots.patient.ratio.sd <- data.frame(row.names = row.names(same_cell.CTRL.boots.patient.ratio.df), SD = rep(0, nrow(same_cell.CTRL.boots.patient.ratio.df)))
for(thePatient in row.names(same_cell.CTRL.boots.patient.ratio.df)){
  same_cell.CTRL.boots.patient.ratio.sd[thePatient,"SD"] <- sd(same_cell.CTRL.boots.patient.ratio.df[thePatient,])
}

# Plot both chrY and CTRL sets together
ggbarplot(m.patient.ratio.bs.median, x = "Patient", y = "value", title = "chrY genes ratio of cells with 0 reads on Y per (male) patient", xlab = "Patient", ylab = "Zero Cell Ratio", fill = "gene.set", palette = "Paired", position = position_dodge()) +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = mean(m.patient.ratio.bs.median[m.patient.ratio.bs.median$gene.set == "CTRL", "value"]), color = brewer.pal(3, "Paired")[2]) +
  geom_hline(yintercept = mean(m.patient.ratio.bs.median[m.patient.ratio.bs.median$gene.set == "chrY", "value"]), color = brewer.pal(3, "Paired")[1]) +
  annotate("text", x = -Inf, y = Inf, label = paste("One-sided two-sample paired rank sum test p-value: ", round(Y.vs.CTRL.pval.bs.median, digits = 5), sep = ""), hjust = -0.1, vjust = 1)

ggsave("mosaicY_macrophages/genes_ratio_per_patient_bootstrapped.same_cells.barplot.pdf")

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
pdf("mosaicY_macrophages/bootstrap_fit_cutoff.same_cells.scatter.pdf")
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
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = curve.fit.cutoff, color = "gold", linewidth = 0.75)
ggsave("mosaicY_macrophages/bootstrap_directionality.same_cells.barplot.pdf")

# Construct ratio diff df
ordered.patient.ratio.bs.median            <- m.egY.patient.ratio
ordered.patient.ratio.bs.median$diff       <- m.egY.patient.ratio$loss.ratio - same_cell.CTRL.boots.patient.ratio.median$value
ordered.patient.ratio.bs.median$loss.ratio <- NULL
ordered.patient.ratio.bs.median            <- merge(ordered.patient.ratio.bs.median, same_cell.CTRL.boots.patient.ratio.sd, by.y = 0, by.x = "Patient")
ordered.patient.ratio.bs.median            <- ordered.patient.ratio.bs.median[order(ordered.patient.ratio.bs.median$diff, decreasing = F),]

# Add in robustness of overrepresentation
m.over_cutoff.robust.boots <- m.robust.boots.df[m.robust.boots.df$Direction == "over",]
m.over_cutoff.robust.boots <- m.over_cutoff.robust.boots[m.over_cutoff.robust.boots[, "Percentage"] >= curve.fit.cutoff, "Patient"]
ordered.patient.ratio.bs.median$Robust <- "No"
ordered.patient.ratio.bs.median[ordered.patient.ratio.bs.median$Patient %in% m.over_cutoff.robust.boots, "Robust"] <- "Yes"

# Plot ratio diff with robustness cutoff
ggbarplot(ordered.patient.ratio.bs.median, x = "Patient", y = "diff", fill = "Robust", title = "Ratio differential of chrY vs. CTRL genes", xlab = "Patient", ylab = 'Difference in "zero cell" ratio (chrY - CTRL)', palette = "Set1") +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_errorbar(aes(ymin=diff-SD, ymax=diff+SD), width=.2,
                position=position_dodge(.9))
ggsave("mosaicY_macrophages/genes_ratio_diff_per_patient_bootstrapped.same_cells.barplot.pdf")



#======================================================================
## Bootstrap the analysis of the CTRL genes in ALL cells for better fidelity
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
  tmp.egNOTY.expr <- as.data.frame(GetAssayData(final.pop.call.from_full.integrated.mac.seurat, slot = "data", assay = "RNA")[CTRL.boots.genes[[n]],])
  
  # Take the average expression of every cell
  tmp.egNOTY.expr <- apply(tmp.egNOTY.expr, 2, mean)
  
  ## Aggregate it back to patient level
  # Get cell names
  egNOTY.cells <- names(tmp.egNOTY.expr)
  
  # Get 10X patient names
  egNOTY.cells.10X    <- egNOTY.cells[grep("plaque", egNOTY.cells)]
  egNOTY.patients.10X <- paste0("P", substr(egNOTY.cells.10X, nchar(egNOTY.cells.10X), nchar(egNOTY.cells.10X)), sep = "")
  
  # Get CEL-seq patient names
  egNOTY.cells.CS    <- egNOTY.cells[grep("plaque", egNOTY.cells, invert = T)]
  egNOTY.patients.CS <- paste0("P", substr(egNOTY.cells.CS, 1, 4), sep = "")
  
  # Merge all patient names
  egNOTY.patients <- c(egNOTY.patients.10X, egNOTY.patients.CS)
  
  # Melt into a frame with patient info
  tmp.m.expr         <- melt(tmp.egNOTY.expr)
  tmp.m.expr$Patient <- egNOTY.patients
  
  # Delete the females
  CTRL.boots.expr.mean.cell[[n]] <- tmp.m.expr[tmp.m.expr$Patient %in% male.patients,]
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
  for(thePatient in unique(male.patients)){
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
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = mean(m.patient.ratio.bs.median[m.patient.ratio.bs.median$gene.set == "CTRL", "value"]), color = brewer.pal(3, "Paired")[2]) +
  geom_hline(yintercept = mean(m.patient.ratio.bs.median[m.patient.ratio.bs.median$gene.set == "chrY", "value"]), color = brewer.pal(3, "Paired")[1]) +
  annotate("text", x = -Inf, y = Inf, label = paste("One-sided two-sample paired rank sum test p-value: ", round(Y.vs.CTRL.pval.bs.median, digits = 5), sep = ""), hjust = -0.1, vjust = 1)

ggsave("mosaicY_macrophages/genes_ratio_per_patient_bootstrapped.barplot.pdf")

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
pdf("mosaicY_macrophages/bootstrap_fit_cutoff.scatter.pdf")
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
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = curve.fit.cutoff, color = "gold", linewidth = 0.75)
ggsave("mosaicY_macrophages/bootstrap_directionality.barplot.pdf")

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
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_errorbar(aes(ymin=diff-SD, ymax=diff+SD), width=.2,
                position=position_dodge(.9))
ggsave("mosaicY_macrophages/genes_ratio_diff_per_patient_bootstrapped.barplot.pdf")



#======================================================================
## Feed mosaic Y metadata back to the Seurat Object
## Construct a metadata column from all mosaic LOY cells (mLOY)
## Only take cells from patients that show a robust difference with controls
# Grab robust patients
robust.mLOY.patients <- ordered.patient.ratio.bs.median[ordered.patient.ratio.bs.median$Robust == "Yes","Patient"]

# Grab zero Y cells form these patients
robust.mLOY.cells <- unlist(egY.zero.cells[names(egY.zero.cells) %in% robust.mLOY.patients])

# Take all cell names and populate a vector of the same length to hold the mLOY metadata
all.cells <- row.names(final.pop.call.integrated.mye.seurat@meta.data)
mLOY.meta <- rep("No", length(all.cells))

# Flip the mLOY positive cells
mLOY.meta[all.cells %in% robust.mLOY.cells] <- "Yes"

# Some stats for sanity checking
length(mLOY.meta[mLOY.meta == "Yes"])
length(mLOY.meta[mLOY.meta == "No"])
length(mLOY.meta)
length(mLOY.meta[mLOY.meta == "Yes"]) / length(mLOY.meta)

# Add the metadata
final.pop.call.integrated.mye.seurat <- AddMetaData(final.pop.call.integrated.mye.seurat, metadata = mLOY.meta, col.name = "mLOY.cell")

# Plot!
customUMAP(object = final.pop.call.integrated.mye.seurat, group.by = "mLOY.cell", pt.size = 1.5, title = "mLOY cells", font.size = 16, legend.pos = "bottom", seed = 666, file.name = "mosaicY_macrophages/mLOY UMAP.pdf", cols = c("grey", "orangered"))
customUMAP(object = final.pop.call.integrated.mye.seurat, pt.size = 1.5, title = "Populations", legend.pos = "none", font.size = 16, seed = 666, file.name = "mosaicY_macrophages/population UMAP.pdf", cols = M.int_refined.pop.colors)

#======================================================================
## Look at proportions of mLOY per population
# Extract the relevant info from the metadata
final.pop.call.integrated.mye.seurat <- AddMetaData(final.pop.call.integrated.mye.seurat, metadata = Idents(final.pop.call.integrated.mye.seurat), col.name = "pop.name")
mLOY.stats <- data.frame( pop.name = final.pop.call.integrated.mye.seurat@meta.data$pop.name, mLOY.cell = final.pop.call.integrated.mye.seurat@meta.data$mLOY.cell, row.names = row.names(final.pop.call.integrated.mye.seurat@meta.data))

# Get the absolute counts
mLOY.per_pop <- vector()
for(thePop in unique(mLOY.stats$pop.name)){
  cat(thePop, "\n")
  
  mLOY.per_pop[thePop] <- sum(mLOY.stats[mLOY.stats$pop.name == thePop, "mLOY.cell"] == "Yes")
}

pop.count <- vector()
for(thePop in unique(mLOY.stats$pop.name)){
  cat(thePop, "\n")
  
  pop.count[thePop] <- nrow(mLOY.stats[mLOY.stats$pop.name == thePop, ])
}

# Get the ratio
mLOY.per_pop.ratio <- vector()
for(thePop in unique(mLOY.stats$pop.name)){
  cat(thePop, "\n")
  
  mLOY.per_pop.ratio[thePop] <- sum(mLOY.stats[mLOY.stats$pop.name == thePop, "mLOY.cell"] == "Yes") / pop.count[thePop]
}

# Concat all info into a df
mLOY.per_pop.df   <- data.frame(total.cells = pop.count, mLOY.cells = mLOY.per_pop, mLOY.ratio = mLOY.per_pop.ratio, pop = names(pop.count))

m.mLOY.per_pop.df  <- melt(mLOY.per_pop.df)
colnames(m.mLOY.per_pop.df) <- c("pop", "count", "cells")
m.mLOY.per_pop.abs <- m.mLOY.per_pop.df[m.mLOY.per_pop.df$count != "mLOY.ratio",]
m.mLOY.per_pop.abs <- m.mLOY.per_pop.abs[c(order(m.mLOY.per_pop.abs[m.mLOY.per_pop.abs$count == "total.cells",]$cells, decreasing = T), 10:18), ]
m.mLOY.per_pop.rel <- m.mLOY.per_pop.df[m.mLOY.per_pop.df$count == "mLOY.ratio", c("pop", "cells")]
m.mLOY.per_pop.rel <- m.mLOY.per_pop.rel[order(m.mLOY.per_pop.rel$cells, decreasing = T), ]
m.mLOY.per_pop.rel$cells <- m.mLOY.per_pop.rel$cells * 100

# Check significance
mLOY.per_pop.pval <- wilcox.test(mLOY.per_pop.df$mLOY.ratio)$p.value

# Plot absolute enrichment
ggbarplot(m.mLOY.per_pop.abs, x = "pop", y = "cells", fill = "count", title = "mLOY cells per population", xlab = "Population", ylab = '# of cells', palette = "Paired") +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1.5/1) +
  annotate("text", x = -Inf, y = Inf, label = paste("Two-sided one-sample rank sum test p-value: ", round(mLOY.per_pop.pval, digits = 3), sep = ""), hjust = -0.1, vjust = 1)

ggsave("mosaicY_macrophages/mLOY.per_pop.bar.pdf")

# Plot relative enrichment
ggbarplot(m.mLOY.per_pop.rel, x = "pop", y = "cells", title = "mLOY cells per population", xlab = "Population", ylab = '% of cells', fill="grey") +
  theme_pubr(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 1.5/1) +
  annotate("text", x = -Inf, y = Inf, label = paste("Two-sided one-sample rank sum test p-value: ", round(mLOY.per_pop.pval, digits = 3), sep = ""), hjust = -0.1, vjust = 1)
ggsave("mosaicY_macrophages/mLOY.per_pop.ratio.bar.pdf")


#======================================================================
## Which Y genes are correlated with foam cell genes, LAM genes, res genes, and inf genes
# Check which Y genes are expressed to begin with
DoHeatmap(final.pop.call.integrated.mye.seurat, features = egY, raster = F, label = F, group.colors = M.int_refined.pop.colors)
ggsave("mosaicY_macrophages/egY.heatmap.pdf")

VlnPlot(final.pop.call.integrated.mye.seurat, features = egY, stack = T, flip = T)
ggsave("mosaicY_macrophages/egY.violin.pdf")

# Remove mac low/unexpressed Y genes
mac_lo.egY <- c("PRORY", "SPRY3", "RBMY2FP", "NLGN4Y", "PPP2R3B", "BCORP1", "GYG2P1", "CRLF2", "TMSB4Y")
mac.egY <- egY[!egY %in% mac_lo.egY]

# Correlation with TREM1
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "TREM1", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.TREM1_cor.pdf",
         width       = 20, 
         height      = 8
)

# Correlation with PLIN2
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "PLIN2", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.PLIN2_cor.pdf",
         width       = 20, 
         height      = 8
)

# Correlation with TREM2
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "TREM2", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.TREM2_cor.pdf",
         width       = 20, 
         height      = 8
)

# Correlation with CD9
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "CD9", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.CD9_cor.pdf",
         width       = 20, 
         height      = 8
)

# Correlation with PLTP
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "PLTP", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.PLTP_cor.pdf",
         width       = 20, 
         height      = 8
)

# Correlation with FOLR2
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "FOLR2", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.FOLR2_cor.pdf",
         width       = 20, 
         height      = 8
)

# Correlation with S100A8
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "S100A8", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.S100A8_cor.pdf",
         width       = 20, 
         height      = 8
)

# Correlation with FCGR3A
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "FCGR3A", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.FCGR3A_cor.pdf",
         width       = 20, 
         height      = 8
)

# Correlation with MX1
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = mac.egY,
         cor.feature = "MX1", 
         ncol        = 7, 
         file.name   = "mosaicY_macrophages/mac_exp.egY.MX1_cor.pdf",
         width       = 20, 
         height      = 8
)


#======================================================================
## Mac expresed Y genes correlation as a gene set

# Retrieve the expression matrix from the seurat object
exprMat <- GetAssayData(final.pop.call.integrated.mye.seurat, "counts")

# Calculate enrichment scores
cells_AUC <- AUCell_run(exprMat, list(egY_expr = mac.egY), BiocParallel::MulticoreParam(5))

# Optimize tresholds
set.seed(333)
par(mfrow=c(1,1)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

cells_assignment$MCSF$aucThr$selected
cells_assignment$GMCSF$aucThr$selected
length(cells_assignment$MCSF$assignment)
length(cells_assignment$GMCSF$assignment)

# Turns out there is no such enrichment


#======================================================================
## Check for correlations with traits
## Glue code needed can be found in final mac populaitons.R

## Process metadata
# Read the file
meta.data <- read.xlsx(file = "raw_data/2021-09-16 AtheroExpress Database Marie.xlsx", sheetIndex = 1, header = T, as.data.frame = T)
dim(meta.data)

# Fix NA values
# Work-around weird 'charToDate()' error when checking if a string is NA while there is a date like structure in the string by forcing those columns to 'factor'
for(theCol in colnames(meta.data)[grep("date", colnames(meta.data), ignore.case = T)]){
  meta.data[, theCol] <- factor(meta.data[, theCol])
}
meta.data[meta.data == "NA"] <- NA

# Filter on quantity of NAs (keep only columns with less than 10% NA)
max.na    <- floor(nrow(meta.data) * 0.10)
meta.data <- meta.data[,colSums(apply(meta.data, 1:2, is.na)) <= max.na]
dim(meta.data)

# Filter out columns with the same value for all rows
meta.data <- meta.data[,apply(meta.data, 2, function(x)length(unique(x))) > 1]
dim(meta.data)

# Filter out columns with the same value or NA for all rows
meta.data <- meta.data[,!(apply(meta.data, 2, function(x)length(unique(x))) == 2 & colSums(apply(meta.data, 1:2, is.na)) > 0)]
dim(meta.data)

## Add metadata to 43p seurat object (as we don't have htis info for the 10X guys)
final.pop.call.integrated.mye.43p.seurat <- subset(final.pop.call.integrated.mye.seurat, subset = Method == "CEL-seq")

# Fetch patient information per cell
md.df <- data.frame(Patient=final.pop.call.integrated.mye.43p.seurat$Patient)

# Expand metadata to all cells
md.df <- merge(md.df, meta.data, by.x = 1, by.y = 2, all.x = T)
md.df$Patient <- NULL

# And add to the Seurat object
for (i in names(md.df)){
  final.pop.call.integrated.mye.43p.seurat <- AddMetaData(final.pop.call.integrated.mye.43p.seurat, md.df[,i], col.name = i)
}

## Let's bin the traits with more than 10 options
binned.meta.data <- meta.data[,1:2]
for(theTrait in colnames(meta.data[,c(-1,-2)])){
  # Keep only traits with over 10 levels
  if(length(levels(factor(meta.data[,theTrait]))) > 10){
    # Keep only numeric traits
    if(check.numeric(meta.data[1,theTrait])){
      cat(paste("Binning:", theTrait,"\n"))
      
      # Define the quartiles of this trait
      q <- quantile(as.numeric(as.character(meta.data[, theTrait])), na.rm = T)
      q.levels <-             paste("n <", round(q[2], digits = 0))
      q.levels <- c(q.levels, paste(round(q[2], digits = 0), "< n <", round(q[3], digits = 0)))
      q.levels <- c(q.levels, paste(round(q[3], digits = 0), "< n <", round(q[4], digits = 0)))
      q.levels <- c(q.levels, paste("n >", round(q[4], digits = 0)))
      
      # Keep a working copy of the trait were we replace NA values
      tmp.df                <- meta.data[,theTrait, drop = F]
      tmp.df[is.na(tmp.df)] <- 0
      
      # Build a results data frame, case by case
      results.df            <- tmp.df
      results.df[,theTrait] <- as.character(results.df[,theTrait])
      
      results.df[as.numeric(as.character(tmp.df[,theTrait]))  < q[2]                                                     , theTrait] <- q.levels[1]
      results.df[as.numeric(as.character(tmp.df[,theTrait])) >= q[2] & as.numeric(as.character(tmp.df[,theTrait])) < q[3], theTrait] <- q.levels[2]
      results.df[as.numeric(as.character(tmp.df[,theTrait])) >= q[3] & as.numeric(as.character(tmp.df[,theTrait])) < q[4], theTrait] <- q.levels[3]
      results.df[as.numeric(as.character(tmp.df[,theTrait])) >= q[4]                                                     , theTrait] <- q.levels[4]
      
      results.df[,theTrait] <- factor(results.df[,theTrait], levels = q.levels)
      
      # Add results back to the main df
      binned.meta.data <- cbind(binned.meta.data, results.df)
    }
  }
}

## Add binned metadata to the seurat object
colnames(binned.meta.data) <- paste0(colnames(binned.meta.data), ".binned", sep = "")
head(binned.meta.data)

# Fetch patient information per cell
md.df <- data.frame(Patient=final.pop.call.integrated.mye.43p.seurat$Patient)

# Expand metadata to all cells
md.df <- merge(md.df, binned.meta.data, by.x = 1, by.y = 2, all.x = T)
md.df$Patient <- NULL

# And add to the Seurat object
for (i in names(md.df)){
  final.pop.call.integrated.mye.43p.seurat <- AddMetaData(final.pop.call.integrated.mye.43p.seurat, md.df[,i], col.name = i)
}


## Loop over the genes and traits
ToI <- c("Symptoms.5G", "AsymptSympt", "Med.statin", "Phenotype", "Age.binned")
male.seurat <- subset(final.pop.call.integrated.mye.43p.seurat, subset = Sex == "male")

for (theGene in mac.egY){
  skip <- "0"
  for(theTrait in ToI){
    cat(theGene, theTrait, "\n")
    nona.seurat <- subset(male.seurat, cells = row.names(male.seurat@meta.data[!is.na(male.seurat@meta.data[,theTrait]),]))
    
    if(which(ToI == theTrait) == 1){
      if(sum(AverageExpression(nona.seurat, features = theGene)$RNA) == 0){
        cat("Skipping",theGene, "...\n\n")
        skip <- "1"
      }else{
        dir.create(paste("mosaicY_macrophages/correlations/", theGene, " Individual trait correlation/", sep = ""), showWarnings = F, recursive = T)
        stratifyByExpression(object = nona.seurat, strat.by = theGene, return.object = F, do.plot = T, onlyUMAP = T, file.name = paste("mosaicY_macrophages/correlations/", theGene, " Individual trait correlation/", theGene, " stratified feature plot.pdf", sep = ""))
      }
    }
    
    if(skip == "0"){
      pv <- kruskal.test(GetAssayData(nona.seurat)[theGene,]~nona.seurat@meta.data[,theTrait])$p.value
      df <- data.frame(x = GetAssayData(nona.seurat)[theGene,], y = nona.seurat@meta.data[,theTrait])
      ggplot(df) + geom_boxplot(aes(y, x, fill = y)) + 
        theme_pubclean(base_size = 32) +
        theme(legend.position = "none") + scale_fill_futurama() +
        ylab(paste(theGene, "Expression", sep = " ")) +
        xlab(theTrait) + 
        annotate(geom = 'text', label = paste("p=",round(pv, digits = 2), sep =""), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = 10, fontface = "bold") + 
        ggtitle(paste(theGene,"correlation with", theTrait, sep = " "))
      dir.create(paste("mosaicY_macrophages/correlations/", theGene, " Individual trait correlation/", sep = ""), showWarnings = F, recursive = T)
      ggsave(filename = paste("mosaicY_macrophages/correlations/", theGene, " Individual trait correlation/", theGene, " correlation with ", theTrait,".pdf", sep = ""))
    }
  }
}


#======================================================================
## PLot some candidates a bit more
mLOY.candidate.culprits <- c("SLC25A6", "CD99", "DDX3Y")
bunchOfCustomPlots(object = final.pop.call.integrated.mye.seurat, 
                   features = mLOY.candidate.culprits, 
                   Vln.draw.names = F, 
                   feature.pt.size = 1.5, 
                   Vln.pt.size = 0.5, 
                   dot.scale = 15, 
                   name = "mosaicY_macrophages/mLOY canditates", 
                   Dot.width = 15,
                   Vln.color = M.int_refined.pop.colors, Vln.width = 20, Vln.height = 6)


#======================================================================
## Check if mLOY status correlates with age
mLOY.age     <- data.frame(mLOY.cell = male.seurat@meta.data$mLOY.cell, Age = male.seurat@meta.data$Age)
not_mLOY.age <- mLOY.age[grep("No", mLOY.age$mLOY.cell), "Age"]
mLOY.age     <- mLOY.age[grep("Yes", mLOY.age$mLOY.cell), "Age"]

mLOY.age     <- sort(mLOY.age)
not_mLOY.age <- sort(not_mLOY.age)

mLOY.age.no_mean <- t.test(Age~mLOY.cell, data = mLOY.age)$estimate[1]
mLOY.age.yes_mean <- t.test(Age~mLOY.cell, data = mLOY.age)$estimate[2]
mLOY.age.pval <- t.test(Age~mLOY.cell, data = mLOY.age)$p.value
ggplot(mLOY.age, aes(x=Age, fill = mLOY.cell)) + 
  geom_histogram(alpha = .5, bins = 10, position = "identity") +
  annotate(geom = 'text', label = paste("p=",round(mLOY.age.pval, digits = 2), sep =""), x = -Inf, y = Inf, hjust = 0, vjust = 1, size = 6, fontface = "bold") +
  theme_pubclean(base_size = 32) + ylab("# of cells") +
  theme(legend.position = "right") + scale_fill_futurama() + ggtitle("Age ~ mLOY") +
  geom_vline(xintercept = c(mLOY.age.no_mean,mLOY.age.yes_mean), color = pal_futurama("planetexpress")(2))
ggsave("mosaicY_macrophages/correlations/mLOY status versus Age.pdf")


