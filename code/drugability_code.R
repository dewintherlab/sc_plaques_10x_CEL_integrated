##==========================================================================================================
## Check drugabilty of target genes (markers + LRI targets) using DGIdb
##=========================================================================================================
# Make directory for drugability results
dir.create("drugability_results/markers_plus_LRI_pairs", showWarnings = FALSE, recursive = T)


## Assemble Myeloid population's target gene lists
# Start with all marker genes
Mye.markers <- list()
for(i in levels(Idents(integrated.mye.seurat))){
  Mye.markers[[i]] <- subset(integrated.mye.seurat.markers, cluster == i & p_val_adj < 0.1)[,"gene"]
  cat(c(i, length(Mye.markers[[i]]), "\n"))
  
}

# Add in cellphone LRI partners (only one of the pair has to be DEG, so we can get some extra hits this way)
for(theIdent in levels(Idents(integrated.mye.seurat))){
  # Convert to the cpdb compatible ident name
  cpdb.ident <- cell_type_map[cell_type_map$seurat.ident == theIdent, "cpdb.ident"]
  
  # Extract ligands and receptor, called in the correct ident only
  cpdb.genes <- levels(unique(c(subset(exp.m.s.cpdb.relevant, ligand.cells   == cpdb.ident)$ligand, 
                         subset(exp.m.s.cpdb.relevant, receptor.cells == cpdb.ident)$receptor)))
  
  # Remove complexes as we will not have them in the DGIdb anyway
  cpdb.genes <- cpdb.genes[grep("complex", cpdb.genes, invert = T)]
  
  # Add to the list
  Mye.markers[[theIdent]] <- unique(c(Mye.markers[[theIdent]], cpdb.genes))
  cat(c(theIdent, length(Mye.markers[[theIdent]]), "\n"))
}

#===========================================================================
#===========================================================================
# Run DGIdb queries
interactionTypes()
Mye.markers.dgi <- list()
for( i in names(Mye.markers)){
  Mye.markers.dgi[[i]] <- queryDGIdb(Mye.markers[[i]])
}
# Save the object
saveRDS(Mye.markers.dgi, "Seurat_Objects/mye_markers_plus_LRI.DGI_object.RDS")
#Mye.markers.dgi <- readRDS("Seurat_Objects/mye_markers_plus_LRI.DGI_object.RDS")

# Write results to disk
for( i in names(Mye.markers)){
  x <- detailedResults(Mye.markers.dgi[[i]])
  x <- merge(x, byGene(Mye.markers.dgi[[i]]))
  x <- merge(x, resultSummary(Mye.markers.dgi[[i]])) 
  x <- x[,c("Gene", "Drug", "InteractionType", "PMIDs", "GeneName", "DistinctDrugCount", "DruggableGeneCategories", "Score")]
  
  x$Score <- as.numeric(x$Score)
  x <- x[order(x$Score, decreasing = T),]
  x <- x[x$DistinctDrugCount > 0 & x$Score > 1,]
  write.table(x, file = paste("drugability_results/markers_plus_LRI_pairs/Cluster ", i, " drug target full table.txt", sep = ""), quote = F, sep = "\t", row.names = F)
}

# Plot number of interactions
for( i in names(Mye.markers.dgi)){
  m <- subset(byGene(Mye.markers.dgi[[i]]),DistinctDrugCount > 0)
  m <- m[order(m$Gene),]
  d <- duplicated(m$Gene)
  m <- m[!d,]
  m <- m[order(m$DistinctDrugCount, decreasing = T),]
  m$Gene <- factor(m$Gene, levels = m$Gene)
  ggplot(m, aes(x = Gene, y = DistinctDrugCount, fill = DistinctDrugCount)) +
           geom_bar(stat = "identity", position = "dodge") +
           coord_flip() +
           scale_fill_gradient(low = "dodgerblue", high = "coral") +
           ylab("Drugs (#)") +
           ggtitle(paste(i, "target/drug interactions")) +
           theme(panel.background = element_blank(),
                 legend.position = "none",
                 axis.ticks.y     = element_blank(),
                 plot.margin      = unit(c(1,1,1,1), "cm"),
                 axis.text        = element_text(size = 18, face = "plain"), 
                 aspect.ratio     = 3/2
    )
  ggsave(paste("drugability_results/markers_plus_LRI_pairs/Cluster ", i, " target drug interaction counts.pdf", sep = ""), width = 20, height = 30)
}

#===========================================================================
# Filter on score
Mye.markers.dgi.filtered <- list()
for( i in names(Mye.markers.dgi)){
  Mye.markers.dgi.filtered[[i]] <- resultSummary(Mye.markers.dgi[[i]])[resultSummary(Mye.markers.dgi[[i]])$Score > 1 ,c("Gene","Drug","Score")]
}

plot.DGI.results()

#===========================================================================
# Filter unique genes (unique to a population)
all.Mye.markers.dgi.filtered <- sort(unlist(lapply(Mye.markers.dgi.filtered, function(x){x <- x[order(x$Gene),]
                                                                                         d <- duplicated(x$Gene)
                                                                                         x <- x[!d,]
                                                                                         return(x$Gene)})))
d <- table(all.Mye.markers.dgi.filtered) == 1
all.Mye.markers.dgi.filtered <- names(d[d])

unique.Mye.markers <- list()
for( i in names(Mye.markers.dgi.filtered)){
  unique.Mye.markers[[i]] <- Mye.markers.dgi.filtered[[i]][Mye.markers.dgi.filtered[[i]]$Gene %in% all.Mye.markers.dgi.filtered,]
}

# Write results to disk
for( i in names(Mye.markers)){
  x       <- unique.Mye.markers[[i]]
  x$Score <- as.numeric(x$Score)
  x       <- x[order(x$Score, decreasing = T),]
  write.table(x, file = paste("drugability_results/markers_plus_LRI_pairs/Cluster ", i, " drug target unique table.txt", sep = ""), quote = F, sep = "\t", row.names = F)
}

plot.DGI.results(dgi.list = unique.Mye.markers, name = "unique potential drugs")


#===========================================================================
#===========================================================================
## Find drugability of broader types (i.e., foam, resident, inflammatory)
## Consolidate the lists
Mye.type.markers <- list()

# Foamy populations
for(theIdent in names(Mye.markers)[grep("Foam", names(Mye.markers))]){
  Mye.type.markers[["Foamy"]] <- unique(c(Mye.type.markers[["Foamy"]], Mye.markers[[theIdent]]))
}

# Resident populations
for(theIdent in names(Mye.markers)[grep("Resident|ABCG", names(Mye.markers))]){
  Mye.type.markers[["Resident"]] <- unique(c(Mye.type.markers[["Resident"]], Mye.markers[[theIdent]]))
}

# Inflammatory populations
for(theIdent in names(Mye.markers)[grep("SELL", names(Mye.markers))]){
  Mye.type.markers[["Inflammatory"]] <- unique(c(Mye.type.markers[["Inflammatory"]], Mye.markers[[theIdent]]))
}

length(Mye.type.markers[["Foamy"]])
length(Mye.type.markers[["Resident"]])
length(Mye.type.markers[["Inflammatory"]])

# Query the DGIdb
Mye.type.markers.dgi <- list()
for( i in names(Mye.type.markers)){
  Mye.type.markers.dgi[[i]] <- queryDGIdb(Mye.type.markers[[i]])
}
# Save the object
saveRDS(Mye.type.markers.dgi, "Seurat_Objects/mye_type_markers_plus_LRI.DGI_object.RDS")
#Mye.type.markers.dgi <- readRDS("Seurat_Objects/mye_type_markers_plus_LRI.DGI_object.RDS")

# Write results to disk
for( i in names(Mye.type.markers)){
  x <- detailedResults(Mye.type.markers.dgi[[i]])
  x <- merge(x, byGene(Mye.type.markers.dgi[[i]]))
  x <- merge(x, resultSummary(Mye.type.markers.dgi[[i]])) 
  x <- x[,c("Gene", "Drug", "InteractionType", "PMIDs", "GeneName", "DistinctDrugCount", "DruggableGeneCategories", "Score")]
  
  x$Score <- as.numeric(x$Score)
  x <- x[order(x$Score, decreasing = T),]
  x <- x[x$DistinctDrugCount > 0 & x$Score > 1,]
  write.table(x, file = paste("drugability_results/markers_plus_LRI_pairs/", i, " drug target full table.txt", sep = ""), quote = F, sep = "\t", row.names = F)
}

## Plot the types
# set up archetype metadata
archetypes <- as.vector(Idents(integrated.mye.seurat))
archetypes[grep("Foamy", archetypes)] <- "Foamy"
archetypes[grep("Resident|ABCG", archetypes)] <- "Resident"
archetypes[grep("SELL", archetypes)] <- "Inflammatory"

# Sanity check
unique(archetypes)

# Add to the seurat object
integrated.mye.seurat <- AddMetaData(integrated.mye.seurat, archetypes, col.name = "archetype")

# Setup the colours
archetype.colors <- c("Resident" = "#528347", "Inflammatory" = "#87482F", "Foamy" = "#9476A3")

# PLot the UMAP
customUMAP(object = integrated.mye.seurat, group.by = "archetype", cols = archetype.colors,
           label = F, shuffle = T, seed = 666, file.name = "drugability_results/markers_plus_LRI_pairs/archetypes.pdf", legend.pos = "none")


# Plot number of interactions
for( i in names(Mye.type.markers.dgi)){
  m <- subset(byGene(Mye.type.markers.dgi[[i]]),DistinctDrugCount > 0)
  m <- m[order(m$Gene),]
  d <- duplicated(m$Gene)
  m <- m[!d,]
  m <- m[order(m$DistinctDrugCount, decreasing = T),]
  m$Gene <- factor(m$Gene, levels = m$Gene)
  ggplot(m, aes(x = Gene, y = DistinctDrugCount, fill = DistinctDrugCount)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    scale_fill_gradient(low = "dodgerblue", high = "coral") +
    ylab("Drugs (#)") +
    ggtitle(paste(i, "target/drug interactions")) +
    theme(panel.background = element_blank(),
          legend.position = "none",
          axis.ticks.y     = element_blank(),
          plot.margin      = unit(c(1,1,1,1), "cm"),
          axis.text        = element_text(size = 18, face = "plain"), 
          aspect.ratio     = 3/2
    )
  ggsave(paste("drugability_results/markers_plus_LRI_pairs/", i, " target drug interaction counts.pdf", sep = ""), width = 20, height = 30)
}


#===========================================================================
# Filter
Mye.type.markers.dgi.filtered <- list()
for( i in names(Mye.type.markers.dgi)){
  Mye.type.markers.dgi.filtered[[i]] <- resultSummary(Mye.type.markers.dgi[[i]])[resultSummary(Mye.type.markers.dgi[[i]])$Score > 1 ,c("Gene","Drug","Score")]
}

plot.DGI.results(dgi.list = Mye.type.markers.dgi.filtered)

#===========================================================================
# Filter unique genes (unique to a population)
all.Mye.type.markers.dgi.filtered <- sort(unlist(lapply(Mye.type.markers.dgi.filtered, function(x){x <- x[order(x$Gene),]
                                                                                                          d <- duplicated(x$Gene)
                                                                                                          x <- x[!d,]
                                                                                                          return(x$Gene)})))
d <- table(all.Mye.type.markers.dgi.filtered) == 1
all.Mye.type.markers.dgi.filtered <- names(d[d])

unique.Mye.type.markers <- list()
for( i in names(Mye.type.markers.dgi.filtered)){
  unique.Mye.type.markers[[i]] <- Mye.type.markers.dgi.filtered[[i]][Mye.type.markers.dgi.filtered[[i]]$Gene %in% all.Mye.type.markers.dgi.filtered,]
}

# Write results to disk
for( i in names(Mye.type.markers)){
  x       <- unique.Mye.type.markers[[i]]
  x$Score <- as.numeric(x$Score)
  x       <- x[order(x$Score, decreasing = T),]
  write.table(x, file = paste("drugability_results/markers_plus_LRI_pairs/" , i, " drug target unique table.txt", sep = ""), quote = F, sep = "\t", row.names = F)
}

plot.DGI.results(dgi.list = unique.Mye.type.markers, name = "unique potential drugs")


#===========================================================================
## Get pathway enrichments
# Per population
ont.unique.Mye.markers <- list()
for( i in names(unique.Mye.markers)){
  theRes                           <- data.frame("gene" = unique(unique.Mye.markers[[i]]$Gene), "avg_log2FC" = 0)
  if(nrow(theRes) > 4){
    ont.unique.Mye.markers[[i]] <- get_ontology(res         = theRes, 
                                                return.data = T, 
                                                plot.top.n  = 10, 
                                                full_GSEA   = F, 
                                                universe    = mye.patients.seurat, 
                                                name        = i,
                                                outdir      = "drugability_results/markers_plus_LRI_pairs")
    }
}


# Per type group
ont.unique.Mye.type.markers <- list()
for( i in names(unique.Mye.type.markers)){
  theRes                           <- data.frame("gene" = unique(unique.Mye.type.markers[[i]]$Gene), "avg_log2FC" = 0)
  ont.unique.Mye.type.markers[[i]] <- get_ontology(res         = theRes, 
                                                   return.data = T, 
                                                   plot.top.n  = 10, 
                                                   full_GSEA   = F, 
                                                   universe    = mye.patients.seurat, 
                                                   name        = i,
                                                   outdir      = "drugability_results/markers_plus_LRI_pairs")
}

##========================================================
## Add pathway info tho the dgi results
##========================================================
# Set up a gene set collection
gsets <- getGeneSets(library = c("H", "C2","C5"))

## Extract relevant pathways
# Foamy pops
pathways.foamy.ont.unique.Mye.type.markers <- list()
for(thePathway in ont.unique.Mye.type.markers$Foamy$Pathways[ont.unique.Mye.type.markers$Foamy$Pathways$padj < 0.05,"name"]){
  cat(paste("Extracting", thePathway, "...\t"))
  if(thePathway %in% names(gsets)){
    pathways.foamy.ont.unique.Mye.type.markers[[thePathway]] <- geneIds(gsets[thePathway])
    cat(paste("Done!\n"))
  }else{cat(paste("Not found!\n"))}
}

# Inflammatory pops
pathways.inf.ont.unique.Mye.type.markers <- list()
for(thePathway in ont.unique.Mye.type.markers$Inflammatory$Pathways[ont.unique.Mye.type.markers$Inflammatory$Pathways$padj < 0.05,"name"]){
  cat(paste("Extracting", thePathway, "...\t"))
  if(thePathway %in% names(gsets)){
    pathways.inf.ont.unique.Mye.type.markers[[thePathway]] <- geneIds(gsets[thePathway])
    cat(paste("Done!\n"))
  }else{cat(paste("Not found!\n"))}
}

# Resident pops. Take padj < 0.1 here or we get only 2 pathways
pathways.res.ont.unique.Mye.type.markers <- list()
for(thePathway in ont.unique.Mye.type.markers$Resident$Pathways[ont.unique.Mye.type.markers$Resident$Pathways$padj < 0.1,"name"]){
  cat(paste("Extracting", thePathway, "...\t"))
  if(thePathway %in% names(gsets)){
    pathways.res.ont.unique.Mye.type.markers[[thePathway]] <- geneIds(gsets[thePathway])
    cat(paste("Done!\n"))
  }else{cat(paste("Not found!\n"))}
}


## Match drugable type markers with pathways
# Foamy pops
unique.Mye.type.markers$Foamy$Pathways <- ""
for(thePathway in names(pathways.foamy.ont.unique.Mye.type.markers)){
  matched.genes <- unlist(pathways.foamy.ont.unique.Mye.type.markers[[thePathway]])[unlist(pathways.foamy.ont.unique.Mye.type.markers[[thePathway]]) %in% unique.Mye.type.markers$Foamy$Gene]
  unique.Mye.type.markers$Foamy[unique.Mye.type.markers$Foamy$Gene %in% matched.genes,"Pathways"] <- paste(unique.Mye.type.markers$Foamy[unique.Mye.type.markers$Foamy$Gene %in% matched.genes,"Pathways"], thePathway, sep = ",")
}

# Inflammatory pops
unique.Mye.type.markers$Inflammatory$Pathways <- ""
for(thePathway in names(pathways.inf.ont.unique.Mye.type.markers)){
  matched.genes <- unlist(pathways.inf.ont.unique.Mye.type.markers[[thePathway]])[unlist(pathways.inf.ont.unique.Mye.type.markers[[thePathway]]) %in% unique.Mye.type.markers$Inflammatory$Gene]
  unique.Mye.type.markers$Inflammatory[unique.Mye.type.markers$Inflammatory$Gene %in% matched.genes,"Pathways"] <- paste(unique.Mye.type.markers$Inflammatory[unique.Mye.type.markers$Inflammatory$Gene %in% matched.genes,"Pathways"], thePathway, sep = ",")
}

# Resident pops
unique.Mye.type.markers$Resident$Pathways <- ""
for(thePathway in names(pathways.res.ont.unique.Mye.type.markers)){
  matched.genes <- unlist(pathways.res.ont.unique.Mye.type.markers[[thePathway]])[unlist(pathways.res.ont.unique.Mye.type.markers[[thePathway]]) %in% unique.Mye.type.markers$Resident$Gene]
  unique.Mye.type.markers$Resident[unique.Mye.type.markers$Resident$Gene %in% matched.genes,"Pathways"] <- paste(unique.Mye.type.markers$Resident[unique.Mye.type.markers$Resident$Gene %in% matched.genes,"Pathways"], thePathway, sep = ",")
}

# Update results on disk
for( i in names(Mye.type.markers)){
  x       <- unique.Mye.type.markers[[i]]
  x$Score <- as.numeric(x$Score)
  x       <- x[order(x$Score, decreasing = T),]
  write.table(x, file = paste("drugability_results/markers_plus_LRI_pairs/" , i, " drug target unique table.txt", sep = ""), quote = F, sep = "\t", row.names = F)
}

## Match pathways with drugable type markers
# Foamy 
intersected.pathways.foamy.ont.unique.Mye.type.markers <- list()
for(thePathway in names(pathways.foamy.ont.unique.Mye.type.markers)){
  intersected.pathways.foamy.ont.unique.Mye.type.markers[[thePathway]] <- unique.Mye.type.markers$Foamy[unique.Mye.type.markers$Foamy$Gene %in% unlist(pathways.foamy.ont.unique.Mye.type.markers[[thePathway]])[unlist(pathways.foamy.ont.unique.Mye.type.markers[[thePathway]]) %in% unique.Mye.type.markers$Foamy$Gene],]
  row.names(intersected.pathways.foamy.ont.unique.Mye.type.markers[[thePathway]]) <- NULL # Reset row names
}

# Inflammatory 
intersected.pathways.inf.ont.unique.Mye.type.markers <- list()
for(thePathway in names(pathways.inf.ont.unique.Mye.type.markers)){
  intersected.pathways.inf.ont.unique.Mye.type.markers[[thePathway]] <- unique.Mye.type.markers$Inflammatory[unique.Mye.type.markers$Inflammatory$Gene %in% unlist(pathways.inf.ont.unique.Mye.type.markers[[thePathway]])[unlist(pathways.inf.ont.unique.Mye.type.markers[[thePathway]]) %in% unique.Mye.type.markers$Inflammatory$Gene],]
  row.names(intersected.pathways.inf.ont.unique.Mye.type.markers[[thePathway]]) <- NULL # Reset row names
}

# Resident 
intersected.pathways.res.ont.unique.Mye.type.markers <- list()
for(thePathway in names(pathways.res.ont.unique.Mye.type.markers)){
  intersected.pathways.res.ont.unique.Mye.type.markers[[thePathway]] <- unique.Mye.type.markers$Resident[unique.Mye.type.markers$Resident$Gene %in% unlist(pathways.res.ont.unique.Mye.type.markers[[thePathway]])[unlist(pathways.res.ont.unique.Mye.type.markers[[thePathway]]) %in% unique.Mye.type.markers$Resident$Gene],]
  row.names(intersected.pathways.foamy.ont.unique.Mye.type.markers[[thePathway]]) <- NULL # Reset row names
}

#========================================
## Plot pathways with drugable markers
dir.create("drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Foamy", showWarnings = F, recursive = T)
dir.create("drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Inflammatory", showWarnings = F, recursive = T)
dir.create("drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Resident", showWarnings = F, recursive = T)

## First, save the lists to disk
# Foam
for(thePathway in names(intersected.pathways.foamy.ont.unique.Mye.type.markers)){
  write.xlsx(x         = intersected.pathways.foamy.ont.unique.Mye.type.markers[[thePathway]], 
             file      = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Foamy_pathway_DGI_intersection.xlsx", 
             sheetName = thePathway, row.names = F, col.names = T, append = T)
}

# Inflammatory
for(thePathway in names(intersected.pathways.inf.ont.unique.Mye.type.markers)){
  write.xlsx(x         = intersected.pathways.inf.ont.unique.Mye.type.markers[[thePathway]], 
             file      = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Inflammatory_pathway_DGI_intersection.xlsx", 
             sheetName = thePathway, row.names = F, col.names = T, append = T)
}

# Resident
for(thePathway in names(intersected.pathways.res.ont.unique.Mye.type.markers)){
  write.xlsx(x         = intersected.pathways.res.ont.unique.Mye.type.markers[[thePathway]], 
             file      = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Resident_pathway_DGI_intersection.xlsx", 
             sheetName = thePathway, row.names = F, col.names = T, append = T)
}

## Now make the plots
plot.DGI.results(dgi.list = intersected.pathways.foamy.ont.unique.Mye.type.markers, save.dir = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Foamy")
plot.DGI.results(dgi.list = intersected.pathways.inf.ont.unique.Mye.type.markers,   save.dir = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Inflammatory")
plot.DGI.results(dgi.list = intersected.pathways.res.ont.unique.Mye.type.markers,   save.dir = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/Resident")

#========================================
## Compare abundance of genes in pathways and scores to make an alternative 'top hit' list
## Count the frequency of genes in various pathways. i.e., identify 'hub genes'
# Foamy
intersected.pathways.foamy.ont.unique.Mye.type.markers.df       <- shape.DGI.gene.pathway.overlap(x = intersected.pathways.foamy.ont.unique.Mye.type.markers)
intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count <- table(t(intersected.pathways.foamy.ont.unique.Mye.type.markers.df))

# Inflammatory
intersected.pathways.inf.ont.unique.Mye.type.markers.df       <- shape.DGI.gene.pathway.overlap(x = intersected.pathways.inf.ont.unique.Mye.type.markers)
intersected.pathways.inf.ont.unique.Mye.type.markers.df.count <- table(t(intersected.pathways.inf.ont.unique.Mye.type.markers.df))

# Resident
intersected.pathways.res.ont.unique.Mye.type.markers.df       <- shape.DGI.gene.pathway.overlap(x = intersected.pathways.res.ont.unique.Mye.type.markers)
intersected.pathways.res.ont.unique.Mye.type.markers.df.count <- table(t(intersected.pathways.res.ont.unique.Mye.type.markers.df))


## Add back score (take the highest scoring if multiple drugs exist for this gene) and abundance (i.e. number of drugs per gene)
## Calculate aggregate score as (count + score*2 + abundance/2) / 3
## i.e. DGIdb score is most important, number of availiable drugs least
# Foamy
intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance           <- data.frame(intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count)
colnames(intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance) <- c("Gene", "Count")
intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance$Count     <- as.numeric(intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance$Count)
intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance           <- merge(intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance,
                                                                                                   concat.DGI.score_and_abundance(intersected.pathways.foamy.ont.unique.Mye.type.markers),
                                                                                                   by = "Gene")
intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance$Aggregate <- (intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance$Count + (intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance$Score*2) + (intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance$Abundance)/2)/3
intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance           <- intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance[order(intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance$Aggregate, decreasing = T),]

# Inflammatory
intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance           <- data.frame(intersected.pathways.inf.ont.unique.Mye.type.markers.df.count)
colnames(intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance) <- c("Gene", "Count")
intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance$Count     <- as.numeric(intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance$Count)
intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance           <- merge(intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance,
                                                                                                   concat.DGI.score_and_abundance(intersected.pathways.inf.ont.unique.Mye.type.markers),
                                                                                                   by = "Gene")
intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance$Aggregate <- (intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance$Count + (intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance$Score*2) + (intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance$Abundance)/2)/3
intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance           <- intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance[order(intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance$Aggregate, decreasing = T),]


# Resident
intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance           <- data.frame(intersected.pathways.res.ont.unique.Mye.type.markers.df.count)
colnames(intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance) <- c("Gene", "Count")
intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance$Count     <- as.numeric(intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance$Count)
intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance           <- merge(intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance,
                                                                                                   concat.DGI.score_and_abundance(intersected.pathways.res.ont.unique.Mye.type.markers),
                                                                                                   by = "Gene")
intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance$Aggregate <- (intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance$Count + (intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance$Score*2) + (intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance$Abundance)/2)/3
intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance           <- intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance[order(intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance$Aggregate, decreasing = T),]


## Make scatterplots
# Foamy
plot.DGI.summary(x        = intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance, 
                 name     = "Foamy populations", 
                 save.dir = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection")

# Inflammatory
plot.DGI.summary(x        = intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance, 
                 name     = "Inflammatory populations", 
                 save.dir = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection")

# Resident
plot.DGI.summary(x        = intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance, 
                 name     = "Resident populations", 
                 save.dir = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection")


## Save the summary dfs
write.xlsx(x         = intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance, 
           file      = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/DGI_hit_summary.xlsx",
           sheetName = "Foamy", 
           row.names = F)
write.xlsx(x         = intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance, 
           file      = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/DGI_hit_summary.xlsx",
           sheetName = "Inflammatory", 
           row.names = F,
           append    = T)
write.xlsx(x         = intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance, 
           file      = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection/DGI_hit_summary.xlsx",
           sheetName = "Resident", 
           row.names = F,
           append    = T)


## Make plots for the top 10 aggregate score genes (top scoring drug per gene only)
## Prep data
tmp.list <- list()

# Foamy
tmp.list[["Foamy Populations"]]               <- Mye.type.markers.dgi.filtered$Foamy[Mye.type.markers.dgi.filtered$Foamy$Gene %in% intersected.pathways.foamy.ont.unique.Mye.type.markers.df.count_score_abundance[1:10,"Gene"], ]
tmp.list[["Foamy Populations"]]               <- tmp.list[["Foamy Populations"]][order(tmp.list[["Foamy Populations"]]$Gene, tmp.list[["Foamy Populations"]]$Score, decreasing = T),]
d                                             <- duplicated(tmp.list[["Foamy Populations"]]$Gene)
tmp.list[["Foamy Populations"]]               <- tmp.list[["Foamy Populations"]][!d,]
tmp.list[["Foamy Populations"]]               <- tmp.list[["Foamy Populations"]][1:10,]
tmp.list[["Foamy Populations"]]$Score         <- as.numeric(tmp.list[["Foamy Populations"]]$Score)
tmp.list[["Foamy Populations"]]               <- tmp.list[["Foamy Populations"]][order(tmp.list[["Foamy Populations"]]$Score, decreasing = T),]

# Inflammatory
tmp.list[["Inflammatory Populations"]]        <- Mye.type.markers.dgi.filtered$Inflammatory[Mye.type.markers.dgi.filtered$Inflammatory$Gene %in% intersected.pathways.inf.ont.unique.Mye.type.markers.df.count_score_abundance[1:10,"Gene"], ]
tmp.list[["Inflammatory Populations"]]        <- tmp.list[["Inflammatory Populations"]][order(tmp.list[["Inflammatory Populations"]]$Gene, tmp.list[["Inflammatory Populations"]]$Score, decreasing = T),]
d                                             <- duplicated(tmp.list[["Inflammatory Populations"]]$Gene)
tmp.list[["Inflammatory Populations"]]        <- tmp.list[["Inflammatory Populations"]][!d,]
tmp.list[["Inflammatory Populations"]]        <- tmp.list[["Inflammatory Populations"]][1:10,]
tmp.list[["Inflammatory Populations"]]$Score  <- as.numeric(tmp.list[["Inflammatory Populations"]]$Score)
tmp.list[["Inflammatory Populations"]]        <- tmp.list[["Inflammatory Populations"]][order(tmp.list[["Inflammatory Populations"]]$Score, decreasing = T),]

# Resident
tmp.list[["Resident Populations"]]            <- Mye.type.markers.dgi.filtered$Resident[Mye.type.markers.dgi.filtered$Resident$Gene %in% intersected.pathways.res.ont.unique.Mye.type.markers.df.count_score_abundance[1:10,"Gene"], ]
tmp.list[["Resident Populations"]]            <- tmp.list[["Resident Populations"]][order(tmp.list[["Resident Populations"]]$Gene, tmp.list[["Resident Populations"]]$Score, decreasing = T),]
d                                             <- duplicated(tmp.list[["Resident Populations"]]$Gene)
tmp.list[["Resident Populations"]]            <- tmp.list[["Resident Populations"]][!d,]
tmp.list[["Resident Populations"]]            <- tmp.list[["Resident Populations"]][1:10,]
tmp.list[["Resident Populations"]]$Score      <- as.numeric(tmp.list[["Resident Populations"]]$Score)
tmp.list[["Resident Populations"]]            <- tmp.list[["Resident Populations"]][order(tmp.list[["Resident Populations"]]$Score, decreasing = T),]

# Make the plots!
plot.DGI.results(dgi.list = tmp.list,
                save.dir  = "drugability_results/markers_plus_LRI_pairs/pathway_drugable_gene_intersection",
                name      = "Aggregate Score top 10 genes")

##========================================================
## Plot the follow up candidates
##========================================================

# All together
bunchOfCustomPlots(object         = integrated.mye.seurat, 
                   features       = c("ITGA4", "RARA", "ALOX5", "AKR1B1"), 
                   assay          = "RNA", 
                   Vln.draw.names = F, 
                   Vln.color      = M.int_refined.pop.colors, 
                   ncol           = 2, 
                   name           = "drugability_results/markers_plus_LRI_pairs/follow_up_targets")

bunchOfCustomPlots(object         = pbmc.monocytes.seurat, 
                   features       = c("ITGA4", "RARA", "ALOX5", "AKR1B1"), 
                   Vln.draw.names = F, 
                   Vln.color      = M_refined.pop.colors, 
                   ncol           = 2, 
                   name           = "drugability_results/markers_plus_LRI_pairs/follow_up_targets_plus_monos")

bunchOfCustomPlots(object         = from_full.integrated.mye.seurat, 
                   features       = c("ITGA4", "RARA", "ALOX5", "AKR1B1"), 
                   Vln.draw.names = F, 
                   Vln.color      = M.int_refined.pop.colors, 
                   ncol           = 2, 
                   name           = "drugability_results/markers_plus_LRI_pairs/follow_up_targets_integrated_monos")


# Per gene stratified
for(theGene in c("ITGA4", "RARA", "ALOX5", "AKR1B1")){
  stratifyByExpression(object        = integrated.mye.seurat, 
                       strat.by      = theGene, 
                       file.name     = paste("drugability_results/markers_plus_LRI_pairs/follow_up_targets ", theGene, ".pdf", sep = ""), 
                       return.object = F, 
                       do.plot       = T, 
                       verbose       = T)
}


## Get pathway information. Search in all pathways, not just enriched ones.
# ITGA4
ITGA4.pathways <- vector()
I <- 1
n <- length(names(gsets))
for(thePathway in names(gsets)){
  cat(paste(I, "of", n, "\t", thePathway, "\n"))
  I <- I + 1
  if("ITGA4" %in% geneIds(gsets[[thePathway]])){ITGA4.pathways <- c(ITGA4.pathways, thePathway)}
}


