# Metadata correlations
dir.create("extra_metadata_binned_noNA", showWarnings = F)

# Read metadata file
md <- read.table(file = "20200204_metadata_selection_18patients_sc_paper.txt", header = T, row.names = 1, sep = "\t", colClasses = "character")
row.names(md) <- substring(row.names(md), first = 3, last = 7)

# Remove empty (NA only) columns
md$crp_avg  <- NULL
md$ep_major <- NULL
md$iph_bin  <- NULL

#========================================================
# Bin the values per column (manually)
colnames(md)

# triglycerides
x <- as.numeric(md$triglyceriden)
y <- x
y[which(x < 1.7)] <- "< 1.7"
y[which(x >= 1.7)] <- ">= 1.7"
md$triglyceriden <- y

# BMI
x <- as.numeric(md$bmi)
y <- x
y[which(x < 19)] <- "< 19"
y[which(x >= 19 & x < 25)] <- ">= 19 & < 25"
y[which(x >= 25 & x < 30)] <- ">= 25 < 30"
y[which(x >= 30)] <- ">= 30"
md$bmi <- y

# total cholesterol
x <- as.numeric(md$totalchol)
y <- x
y[which(x < 5)] <- "< 5"
y[which(x >= 5 & x < 8)] <- ">= 5 & < 8"
y[which(x >= 8)] <- ">= 8"
md$totalchol <- y

# LDL
x <- as.numeric(md$ldl)
y <- x
y[which(x < 3)] <- "< 3"
y[which(x >= 3)] <- ">= 3"
md$ldl <- y

# HDL
x <- as.numeric(as.character(md$hdl))
x2 <- md$Gender
y <- x
y[which(x < 1 & x2 == "Male")] <- "< 1 (Males)"
y[which(x < 1.2 & x2 == "Female")] <- "< 1.2 (Females)"
y[which(x >= 1 & x2 == "Male")] <- ">= 1 (Males)"
y[which(x >= 1.2 & x2 == "Female")] <- ">= 1.2 (Females)"
md$hdl <- y

# GFR
x <- as.numeric(md$gfr_mdrd)
y <- x
y[which(x >= 30 & x < 45)] <- ">=30 & < 45"
y[which(x >= 45 & x < 60)] <- ">= 45 & < 60"
y[which(x >= 60 & x < 90)] <- ">= 60 & < 90"
y[which(x >= 90)] <- ">= 90"
md$gfr_mdrd <- y

# pack years
x <- as.numeric(md$epackyearssmoking)
y <- x
y[which(x <= 1)] <- "<= 1"
y[which(x > 1 & x <= 10)] <- "> 1 & <= 10"
y[which(x > 10)] <- "> 10"
md$epackyearssmoking <- y

# Age
x <- as.numeric(md$age)
y <- x
y[which(x < 60)] <- "< 60"
y[which(x >= 60 & x < 70)] <- ">= 60 & < 70"
y[which(x >= 70 & x < 80)] <- ">= 70 & < 80"
y[which(x >= 80)] <- ">= 80"
md$age <- y

head(md)

#========================================================
#Remove NAs
md[is.na(md)] <- "Unknown"
md[] <- lapply(md, factor)

# Fetch patient information per cell
md.df <- data.frame(Patient=v3.all.seur.combined@meta.data$Patient)

# Expand metadata to all cells
md.df <- merge(md.df, md, by.x = 1, by.y = 0, all.x = T)
md.df$Patient <- NULL

# Add to the seurat object
for (i in names(md.df)){
  v3.all.seur.combined <- AddMetaData(v3.all.seur.combined, md.df[,i], col.name = i)
}

# Visualize
for(i in colnames(md.df)){
  TSNEPlot(v3.all.seur.combined, group.by = i,  pt.size = 2, label.size = 10) + theme(panel.background = element_blank(),
                                                                                      text             = element_text(size = 16, face = "bold"),
                                                                                      axis.text        = element_text(size = 14, face = "plain")
  )
  ggsave(paste("extra_metadata_binned_noNA/", i, "_tSNE.pdf", sep = ""), width = 15, height = 15)
}


#=============================================================================================
# Prepare data frame containing the number of cells per cluster
# Populate with total cells per patient
num_cells.per_patient <- data.frame()
for(i in unique(v3.all.seur.combined$Patient)){
  num_cells.per_patient[i, "Total"] <- sum(v3.all.seur.combined$Patient == i)
}
# Display
num_cells.per_patient 

# Add cells per cluster
for(i in unique(v3.all.seur.combined$Patient)){
  for(j in unique(Idents(v3.all.seur.combined))){
    num_cells.per_patient[i, j] <- sum(v3.all.seur.combined$Patient == i & Idents(v3.all.seur.combined) == j)
  }
}

# Order on patient number
num_cells.per_patient <- num_cells.per_patient[order(row.names(num_cells.per_patient)),]

# Sanity check
cat("Rows are sane: ", sum((rowSums(num_cells.per_patient) - num_cells.per_patient$Total) == num_cells.per_patient$Total) == length(unique(v3.all.seur.combined$Patient)), "\n")
cat("Columns are sane: ", sum(colSums(num_cells.per_patient[,2:dim(num_cells.per_patient)[2]])[order(colnames(num_cells.per_patient[,2:dim(num_cells.per_patient)[2]]))] == cell.numbers[order(names(cell.numbers))]) == length(unique(Idents(v3.all.seur.combined))), "\n")

# Display
num_cells.per_patient

# Remove outlier patient 4440 (only 7 cells total o_O)
num_cells.per_patient <- num_cells.per_patient[which(row.names(num_cells.per_patient) != "4440"),]

# move 'Totals column to it's own frame
total_cells.per_patient <- num_cells.per_patient[,"Total", drop = F]
num_cells.per_patient$Total <- NULL


#=============================================================================================
# Check distribution of metadata vers over the clusters
# Split per factor level
md.vars <- list()
for (theVar in colnames(md)){
  md.vars[[theVar]] <- list()
  for (i in levels(md[,theVar])){
    md.vars[[theVar]][[i]] <- num_cells.per_patient[which(md[,theVar] == i),]
    # Clean up NAs
    md.vars[[theVar]][[i]] <- md.vars[[theVar]][[i]][rowSums(is.na(md.vars[[theVar]][[i]])) != ncol(md.vars[[theVar]][[i]]), ]
  }
}

# Remove Unknowns
for (theVar in colnames(md)){
  md.vars[[theVar]][["Unknown"]] <- NULL
}

# Combine counts per type
summed.md.vars <- list()
t.summed.md.vars <- list()
for (theVar in colnames(md)){
  summed.md.vars[[theVar]]   <- data.frame(row.names = colnames(num_cells.per_patient))
  t.summed.md.vars[[theVar]] <- data.frame(row.names = colnames(num_cells.per_patient))
  for (i in levels(md[,theVar])[grep("Unknown", levels(md[,theVar]), invert = T)]){
    summed.md.vars[[theVar]][,i]  <- apply(md.vars[[theVar]][[i]], 2, sum)
    t.summed.md.vars[[theVar]][,i] <- apply(md.vars[[theVar]][[i]], 2, sum)
  }
  summed.md.vars[[theVar]] <- t(summed.md.vars[[theVar]])
}

# Run Fisher exact tests
fisher.summed.md.vars <- list()
for (theVar in colnames(md)){
  for (i in levels(md[,theVar])[grep("Unknown", levels(md[,theVar]), invert = T)]){ # Loop over all metadata
    for (j in 1:(dim(summed.md.vars[[theVar]])[2]-1)){ #Loop over levels 1:n-1
      sublist <- list()
      for (k in (j+1):dim(summed.md.vars[[theVar]])[2]){ #Loop over levels 2:n
        sublist[[colnames(summed.md.vars[[theVar]])[k]]] <- fisher.test(rbind(summed.md.vars[[theVar]][,j], summed.md.vars[[theVar]][,k]), simulate.p.value = T)
      }
      fisher.summed.md.vars[[theVar]][[colnames(summed.md.vars[[theVar]])[j]]] <- sublist
    }
  }
}

t.fisher.summed.md.vars <- list()
for (theVar in colnames(md)){
  for (i in levels(md[,theVar])[grep("Unknown", levels(md[,theVar]), invert = T)]){ # Loop over all metadata
    for (j in 1:(dim(t.summed.md.vars[[theVar]])[2]-1)){ #Loop over levels 1:n-1
      sublist <- list()
      t.sublist <- list()
      for (k in (j+1):dim(t.summed.md.vars[[theVar]])[2]){ #Loop over levels 2:n
        if(length(unique(t.summed.md.vars[[theVar]][,j])) != 1 & length(unique(t.summed.md.vars[[theVar]][,k])) != 1){
          t.sublist[[colnames(t.summed.md.vars[[theVar]])[k]]] <- fisher.test(rbind(t.summed.md.vars[[theVar]][,j], t.summed.md.vars[[theVar]][,k]), simulate.p.value = T)
        }
      }
      t.fisher.summed.md.vars[[theVar]][[colnames(t.summed.md.vars[[theVar]])[j]]] <- t.sublist
    }
  }
}

# Extract p values
pvals.fisher.summed.md.vars <- list()
for (theVar in colnames(md)){
  for (i in 1:(dim(summed.md.vars[[theVar]])[2]-1)){
    for (j in (i+1):dim(summed.md.vars[[theVar]])[2]){
      pvals.fisher.summed.md.vars[[theVar]] <- rbind(pvals.fisher.summed.md.vars[[theVar]] , 
                                                    data.frame(First   = colnames(num_cells.per_patient)[i],
                                                               Second  = colnames(num_cells.per_patient)[j],
                                                               p.value = fisher.summed.md.vars[[theVar]][[colnames(num_cells.per_patient)[i]]][[colnames(num_cells.per_patient)[j]]]$p.value
                                                    )
      )
    }
  }
}

t.pvals.fisher.summed.md.vars <- list()
for (theVar in colnames(md)){
  for (i in 1:(dim(t.summed.md.vars[[theVar]])[2]-1)){
    for (j in (i+1):dim(t.summed.md.vars[[theVar]])[2]){
      if(length(unique(t.summed.md.vars[[theVar]][,i])) != 1 & length(unique(t.summed.md.vars[[theVar]][,j])) != 1){
        t.pvals.fisher.summed.md.vars[[theVar]] <- rbind(t.pvals.fisher.summed.md.vars[[theVar]] , 
                                                        data.frame(First   = colnames(t.summed.md.vars[[theVar]])[i],
                                                                   Second  = colnames(t.summed.md.vars[[theVar]])[j],
                                                                   p.value = t.fisher.summed.md.vars[[theVar]][[colnames(t.summed.md.vars[[theVar]])[i]]][[colnames(t.summed.md.vars[[theVar]])[j]]]$p.value
                                                        )
        )
      }
    }
  }
}

#Adjust p value for multiple testing
for (theVar in colnames(md)){
  pvals.fisher.summed.md.vars[[theVar]]$padj   <- p.adjust(pvals.fisher.summed.md.vars[[theVar]]$p.value, method = "BH")
  t.pvals.fisher.summed.md.vars[[theVar]]$padj <- p.adjust(t.pvals.fisher.summed.md.vars[[theVar]]$p.value, method = "BH")
}

# Check for significance
lapply(pvals.fisher.summed.md.vars, function(x) x[which(x$padj < 0.05),])
lapply(t.pvals.fisher.summed.md.vars, function(x) x[which(x$padj < 0.05),])

# Build a plot based on percentages to visualize
perc.summed.md.vars <- list()
for (theVar in colnames(md)){
  perc.summed.md.vars[[theVar]] <- prop.table(summed.md.vars[[theVar]], 2)
}

t.perc.summed.md.vars <- list()
for (theVar in colnames(md)){
  t.perc.summed.md.vars[[theVar]] <- prop.table(t(summed.md.vars[[theVar]]), 2)
}

for (theVar in colnames(md)){
  d <- melt(perc.summed.md.vars[[theVar]])
  colnames(d) <- c("Type", "Cluster", "Percentage")
  d$Percentage <- d$Percentage * 100
  
  ggplot(d, aes(x = Cluster, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(type = "qual", palette = "Set3", aesthetics = "fill") +
    ylab("Percentage of Cells") +
    ggtitle(paste("Main clusters - ", theVar, sep = "")) +
    theme_light() +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(size = 12,
                                     angle = 45, 
                                     hjust = 1),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          panel.grid = element_blank(), 
          aspect.ratio = 1/2
    )
  ggsave(paste("extra_metadata_binned_noNA/main_clusters_perc.", theVar, ".pdf"))
}

for (theVar in colnames(md)){
  d <- melt(t.perc.summed.md.vars[[theVar]])
  colnames(d) <- c("Type", "Cluster", "Percentage")
  d$Percentage <- d$Percentage * 100
  
  ggplot(d, aes(x = Cluster, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(type = "qual", palette = "Set3", aesthetics = "fill") +
    ylab("Percentage of Cells") +
    ggtitle(paste("Main clusters - ", theVar, sep = "")) +
    theme_light() +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(size = 12,
                                     angle = 45, 
                                     hjust = 1),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          panel.grid = element_blank(), 
          aspect.ratio = 1/2
    )
  ggsave(paste("extra_metadata_binned_noNA/transposed_main_clusters_perc.", theVar, ".pdf", sep = ""))
}


#===========================================================================
#===========================================================================
#Let's look at some sub-clusters
# CD4, CD8, Endo subsets: not enough cells in a lot of patients!
#Macs reasonably OK so let's check those

# Extract number of cells per patient
# Update our object to seurat v3
v3.all.seur.combined.M_clusters <- UpdateSeuratObject(all.seur.combined.M_clusters)

# Populate with total cells per patient
num_cells.per_patient.M_clusters <- data.frame()
for(i in unique(v3.all.seur.combined.M_clusters$Patient)){
  num_cells.per_patient.M_clusters[i, "Total"] <- sum(v3.all.seur.combined.M_clusters$Patient == i)
}
# Display
num_cells.per_patient.M_clusters

# Add cells per cluster
for(i in unique(v3.all.seur.combined.M_clusters$Patient)){
  for(j in unique(Idents(v3.all.seur.combined.M_clusters))){
    num_cells.per_patient.M_clusters[i, j] <- sum(v3.all.seur.combined.M_clusters$Patient == i & Idents(v3.all.seur.combined.M_clusters) == j)
  }
}

# Order on patient number
num_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[order(row.names(num_cells.per_patient.M_clusters)),]

# Display
num_cells.per_patient.M_clusters

# Remove outliers patient 4440 and 4477
num_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[which(row.names(num_cells.per_patient.M_clusters) != "4440"),]
num_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[which(row.names(num_cells.per_patient.M_clusters) != "4477"),]

# Remove 'T-cell contaminent' subcluster 4
num_cells.per_patient.M_clusters$Total <- num_cells.per_patient.M_clusters$Total - num_cells.per_patient.M_clusters$`4`
num_cells.per_patient.M_clusters$`4` <- NULL

# move 'Totals column to it's own frame
total_cells.per_patient.M_clusters <- num_cells.per_patient.M_clusters[,"Total", drop = F]
num_cells.per_patient.M_clusters$Total <- NULL

#=============================================================================================
# Check distribution of metadata vers over the clusters
# Split per factor level
md.vars.M_clusters <- list()
for (theVar in colnames(md)){
  md.vars.M_clusters[[theVar]] <- list()
  for (i in levels(md[,theVar])){
    md.vars.M_clusters[[theVar]][[i]] <- num_cells.per_patient.M_clusters[which(md[,theVar] == i),]
    # Clean up NAs
    md.vars.M_clusters[[theVar]][[i]] <- md.vars.M_clusters[[theVar]][[i]][rowSums(is.na(md.vars.M_clusters[[theVar]][[i]])) != ncol(md.vars.M_clusters[[theVar]][[i]]), ]
  }
}

# Remove Unknowns
for (theVar in colnames(md)){
  md.vars[[theVar]][["Unknown"]] <- NULL
}

# Combine counts per type
summed.md.vars.M_clusters <- list()
t.summed.md.vars.M_clusters <- list()
for (theVar in colnames(md)){
  summed.md.vars.M_clusters[[theVar]]   <- data.frame(row.names = colnames(num_cells.per_patient.M_clusters))
  t.summed.md.vars.M_clusters[[theVar]] <- data.frame(row.names = colnames(num_cells.per_patient.M_clusters))
  for (i in levels(md[,theVar])[grep("Unknown", levels(md[,theVar]), invert = T)]){
    summed.md.vars.M_clusters[[theVar]][,i]  <- apply(md.vars.M_clusters[[theVar]][[i]], 2, sum)
    t.summed.md.vars.M_clusters[[theVar]][,i] <- apply(md.vars.M_clusters[[theVar]][[i]], 2, sum)
  }
  summed.md.vars.M_clusters[[theVar]] <- t(summed.md.vars.M_clusters[[theVar]])
}

# Run Fisher exact tests
fisher.summed.md.vars.M_clusters <- list()
for (theVar in colnames(md)){
  for (i in levels(md[,theVar])[grep("Unknown", levels(md[,theVar]), invert = T)]){ # Loop over all metadata
    for (j in 1:(dim(summed.md.vars.M_clusters[[theVar]])[2]-1)){ #Loop over levels 1:n-1
      sublist <- list()
      for (k in (j+1):dim(summed.md.vars.M_clusters[[theVar]])[2]){ #Loop over levels 2:n
        sublist[[colnames(summed.md.vars.M_clusters[[theVar]])[k]]] <- fisher.test(summed.md.vars.M_clusters[[theVar]][,j], summed.md.vars.M_clusters[[theVar]][,k])
      }
      fisher.summed.md.vars.M_clusters[[theVar]][[colnames(summed.md.vars.M_clusters[[theVar]])[j]]] <- sublist
    }
  }
}

t.fisher.summed.md.vars.M_clusters <- list()
for (theVar in colnames(md)){
  for (i in levels(md[,theVar])[grep("Unknown", levels(md[,theVar]), invert = T)]){ # Loop over all metadata
    for (j in 1:(dim(t.summed.md.vars.M_clusters[[theVar]])[2]-1)){ #Loop over levels 1:n-1
      sublist <- list()
      t.sublist <- list()
      for (k in (j+1):dim(t.summed.md.vars.M_clusters[[theVar]])[2]){ #Loop over levels 2:n
        if(length(unique(t.summed.md.vars.M_clusters[[theVar]][,j])) != 1 & length(unique(t.summed.md.vars.M_clusters[[theVar]][,k])) != 1){
          t.sublist[[colnames(t.summed.md.vars.M_clusters[[theVar]])[k]]] <- fisher.test(t.summed.md.vars.M_clusters[[theVar]][,j], t.summed.md.vars.M_clusters[[theVar]][,k]) 
        }
      }
      t.fisher.summed.md.vars.M_clusters[[theVar]][[colnames(t.summed.md.vars.M_clusters[[theVar]])[j]]] <- t.sublist
    }
  }
}

# Extract p values
pvals.fisher.summed.md.vars.M_clusters <- list()
for (theVar in colnames(md)){
  for (i in 1:(dim(summed.md.vars.M_clusters[[theVar]])[2]-1)){
    for (j in (i+1):dim(summed.md.vars.M_clusters[[theVar]])[2]){
      pvals.fisher.summed.md.vars.M_clusters[[theVar]] <- rbind(pvals.fisher.summed.md.vars.M_clusters[[theVar]] , 
                                                               data.frame(First   = colnames(num_cells.per_patient.M_clusters)[i],
                                                                          Second  = colnames(num_cells.per_patient.M_clusters)[j],
                                                                          p.value = fisher.summed.md.vars.M_clusters[[theVar]][[colnames(num_cells.per_patient.M_clusters)[i]]][[colnames(num_cells.per_patient.M_clusters)[j]]]$p.value
                                                               )
      )
    }
  }
}

t.pvals.fisher.summed.md.vars.M_clusters <- list()
for (theVar in colnames(md)){
  for (i in 1:(dim(t.summed.md.vars.M_clusters[[theVar]])[2]-1)){
    for (j in (i+1):dim(t.summed.md.vars.M_clusters[[theVar]])[2]){
      if(length(unique(t.summed.md.vars.M_clusters[[theVar]][,i])) != 1 & length(unique(t.summed.md.vars.M_clusters[[theVar]][,j])) != 1){
        t.pvals.fisher.summed.md.vars.M_clusters[[theVar]] <- rbind(t.pvals.fisher.summed.md.vars.M_clusters[[theVar]] , 
                                                                   data.frame(First   = colnames(t.summed.md.vars.M_clusters[[theVar]])[i],
                                                                              Second  = colnames(t.summed.md.vars.M_clusters[[theVar]])[j],
                                                                              p.value = t.fisher.summed.md.vars.M_clusters[[theVar]][[colnames(t.summed.md.vars.M_clusters[[theVar]])[i]]][[colnames(t.summed.md.vars.M_clusters[[theVar]])[j]]]$p.value
                                                                   )
        )
      }
    }
  }
}

#Adjust p value for multiple testing
for (theVar in colnames(md)){
  pvals.fisher.summed.md.vars.M_clusters[[theVar]]$padj   <- p.adjust(pvals.fisher.summed.md.vars.M_clusters[[theVar]]$p.value, method = "BH")
  t.pvals.fisher.summed.md.vars.M_clusters[[theVar]]$padj <- p.adjust(t.pvals.fisher.summed.md.vars.M_clusters[[theVar]]$p.value, method = "BH")
}

# Check for significance
lapply(pvals.fisher.summed.md.vars.M_clusters, function(x) x[which(x$padj < 0.05),])
lapply(pvals.fisher.summed.md.vars.M_clusters, function(x) x[which(x$p.value < 0.05),])
### Nothing significant ###

lapply(t.pvals.fisher.summed.md.vars.M_clusters, function(x) x[which(x$padj < 0.05),])
lapply(t.pvals.fisher.summed.md.vars.M_clusters, function(x) x[which(x$p.value < 0.05),])
### Nothing significant ###

# Build a plot based on percentages to visualize
perc.summed.md.vars.M_clusters <- list()
for (theVar in colnames(md)){
  perc.summed.md.vars.M_clusters[[theVar]] <- prop.table(summed.md.vars.M_clusters[[theVar]], 2)
}

t.perc.summed.md.vars.M_clusters <- list()
for (theVar in colnames(md)){
  t.perc.summed.md.vars.M_clusters[[theVar]] <- prop.table(t(summed.md.vars.M_clusters[[theVar]]), 2)
}

for (theVar in colnames(md)){
  d <- melt(perc.summed.md.vars.M_clusters[[theVar]])
  colnames(d) <- c("Type", "Cluster", "Percentage")
  d$Percentage <- d$Percentage * 100
  
  ggplot(d, aes(x = Cluster, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(type = "qual", palette = "Set3", aesthetics = "fill") +
    ylab("Percentage of Cells") +
    ggtitle(paste("Myeloid clusters - ", theVar, sep = "")) +
    theme_light() +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(size = 12,
                                     angle = 45, 
                                     hjust = 1),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          panel.grid = element_blank(), 
          aspect.ratio = 1/2
    )
  ggsave(paste("extra_metadata_binned_noNA/M_clusters_perc.", theVar, ".pdf"))
}

for (theVar in colnames(md)){
  d <- melt(t.perc.summed.md.vars.M_clusters[[theVar]])
  colnames(d) <- c("Type", "Cluster", "Percentage")
  d$Percentage <- d$Percentage * 100
  d$Cluster <- as.factor(d$Cluster)
  d$Type <- as.factor(d$Type)
  
  ggplot(d, aes(x = Cluster, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_brewer(type = "qual", palette = "Set3", aesthetics = "fill") +
    ylab("Percentage of Cells") +
    ggtitle(paste("Myeloid clusters - ", theVar, sep = "")) +
    theme_light() +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(size = 12,
                                     angle = 45, 
                                     hjust = 1),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          panel.grid = element_blank(), 
          aspect.ratio = 1/2
    )
  ggsave(paste("extra_metadata_binned_noNA/transposed_M_clusters_perc.", theVar, ".pdf", sep =""))
}

