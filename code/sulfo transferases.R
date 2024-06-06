#======================================================================
## Analyze sulfo transferases in the mono + mac cells
#======================================================================
dir.create("sulfo transferase results/shortlist", showWarnings = F, recursive = T)
dir.create("sulfo transferase results/all", showWarnings = F, recursive = T)

## First check the specifically requested short list of genes
sultran.genes.shortlist <- c("CHST2", "CHST15", "SELPLG")

DefaultAssay(final.pop.call.from_full.integrated.mye.seurat) <- "RNA"
customUMAP(object = final.pop.call.from_full.integrated.mye.seurat, pt.size = 1.5, cols = M.int_refined.pop.colors, legend.pos = "right", file.name = "sulfo transferase results/Population UMAP.pdf", plot.width = 15, shuffle = T)
bunchOfCustomPlots(object          = final.pop.call.from_full.integrated.mye.seurat, 
                   features        = sultran.genes.shortlist, 
                   Vln.draw.names  = F, 
                   Vln.pt.size     = 0.5, 
                   feature.pt.size = 1.5, 
                   Vln.width       = 17, 
                   Vln.height      = 5, 
                   name            = "sulfo transferase results/shortlist/SLAN related genes", 
                   Vln.color       = M.int_refined.pop.colors, 
                   dot.scale       = 15, 
                   Dot.width       = 15)


## Now check all sulfo transferases
# Build the list (taken from wikipedia)
sultran.genes <- list()
sultran.genes[["carbohydrate"]]       <- row.names(final.pop.call.from_full.integrated.mye.seurat)[grep("^CHST",   row.names(final.pop.call.from_full.integrated.mye.seurat))]
sultran.genes[["galactose-2-O"]]      <- row.names(final.pop.call.from_full.integrated.mye.seurat)[grep("^GAL3ST", row.names(final.pop.call.from_full.integrated.mye.seurat))]
sultran.genes[["heparansulfate-2-0"]] <- row.names(final.pop.call.from_full.integrated.mye.seurat)[grep("^HS2ST1", row.names(final.pop.call.from_full.integrated.mye.seurat))]
sultran.genes[["heparansulfate-3-0"]] <- row.names(final.pop.call.from_full.integrated.mye.seurat)[grep("^HS3ST",  row.names(final.pop.call.from_full.integrated.mye.seurat))]
sultran.genes[["heparansulfate-6-0"]] <- row.names(final.pop.call.from_full.integrated.mye.seurat)[grep("^HS6ST",  row.names(final.pop.call.from_full.integrated.mye.seurat))]
sultran.genes[["N-deacetylase"]]      <- row.names(final.pop.call.from_full.integrated.mye.seurat)[grep("^NDST",   row.names(final.pop.call.from_full.integrated.mye.seurat))]
sultran.genes[["tyrosylprotein"]]     <- row.names(final.pop.call.from_full.integrated.mye.seurat)[grep("^TPST",   row.names(final.pop.call.from_full.integrated.mye.seurat))]
sultran.genes[["sulfo"]]              <- row.names(final.pop.call.from_full.integrated.mye.seurat)[grep("^SULT",   row.names(final.pop.call.from_full.integrated.mye.seurat))]
sultran.genes

# Plot per category
for(theType in names(sultran.genes)){
  cat(theType, "\n")
  dir.create(paste("sulfo transferase results/all/",theType, " sulfotransferases/", sep =""), showWarnings = F, recursive = T)
  for(theGene in sultran.genes[[theType]]){
    cat("\t", theGene, "...\n")
    tryCatch(expr = bunchOfCustomPlots(object          = final.pop.call.from_full.integrated.mye.seurat, 
                                       features        = theGene, 
                                       Vln.draw.names  = F, 
                                       Vln.pt.size     = 0.5, 
                                       feature.pt.size = 1.5, 
                                       Vln.width       = 10, 
                                       Vln.height      = 5, 
                                       name            = paste("sulfo transferase results/all/",theType, " sulfotransferases/", theGene, sep =""), 
                                       Vln.color       = M.int_refined.pop.colors, 
                                       dot.scale       = 15, 
                                       Dot.width       = 15),error = function(e){cat("\t\tFAILED!\n")}, warning = function(e) { "\t\tWARNING!"}, finally = {cat("\tOK!\n")}
            
    )
  }
}

