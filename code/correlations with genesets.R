#======================================================================
## Correlation of CELC4A with resident and other genes.
## In the paper, they compare per patient expression of these genes, in a bulk transcriptome data set. So withouth cell type info
## First we emulate this approach, then we go more granular and look per cell type
#======================================================================
dir.create("gene_set_correlations", showWarnings = F)

# Load the monaco resident gene set + CLEC4A
monaco.res.set <- c("CD163", "CSF1R", "CD200R1", "CLEC7A", "TLR7", "SCARB1", "CD68", "MSR1", "MRC1", "MS4A4A")

# Check availabilty in our data: 0 == All available
sum(!monaco.res.set %in% row.names(final.pop.call.integrated.mye.seurat))


#======================================================================
## Emulate the paper (so (pseudo-)bulk seq per patient)
# Check Monaco's res genes
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = monaco.res.set,
         cor.feature = "CLEC4A", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/all_cells.res_genes.CLEC4A.pdf",
         width       = 20, 
         height      = 8
         )

# Check Monaco's res genes but correlate with CD14
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = monaco.res.set,
         cor.feature = "CD14", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/all_cells.res_genes.CD14.pdf",
         width       = 20, 
         height      = 8
)

# Check Monaco's res genes but correlate with CD8A
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = monaco.res.set,
         cor.feature = "CD8A", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/all_cells.res_genes.CD8A.pdf",
         width       = 20, 
         height      = 8
)

# Check Monaco's res genes but correlate with TIMP1
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = monaco.res.set,
         cor.feature = "TIMP1", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/all_cells.res_genes.TIMP1.pdf",
         width       = 20, 
         height      = 8
)

# Check Monaco's res genes but correlate with IRF5
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = monaco.res.set,
         cor.feature = "IRF5", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/all_cells.res_genes.IRF5.pdf",
         width       = 20, 
         height      = 8
)

# Do the id genes
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "CLEC4A", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/all_cells.id_genes.CLEC4A.pdf",
         width       = 20, 
         height      = 20
         )

# Do the id genes but correlate them with CD14
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "CD14", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/all_cells.id_genes.CD14.pdf",
         width       = 20, 
         height      = 20
)

# Do the id genes but correlate them with CD8A
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "CD8A", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/all_cells.id_genes.CD8A.pdf",
         width       = 20, 
         height      = 20
)

# Do the id genes but correlate them with TIMP1
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "TIMP1", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/all_cells.id_genes.TIMP1.pdf",
         width       = 20, 
         height      = 20
)

# Do the id genes but correlate them with IRF5
plot.cor(object      = final.pop.call.integrated.full.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "IRF5", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/all_cells.id_genes.IRF5.pdf",
         width       = 20, 
         height      = 20
)

#======================================================================
## Now check only in the macrophages (as a whole)
# Check Monaco's res genes
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = monaco.res.set,
         cor.feature = "CLEC4A", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/macrophages.res_genes.CLEC4A.pdf",
         width       = 20, 
         height      = 8
)

# Check Monaco's res genes but correlate with CD14
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = monaco.res.set,
         cor.feature = "CD14", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/macrophages.res_genes.CD14.pdf",
         width       = 20, 
         height      = 8
)

# Check Monaco's res genes but correlate with CASP1
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = monaco.res.set,
         cor.feature = "CASP1", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/macrophages.res_genes.CASP1.pdf",
         width       = 20, 
         height      = 8
)

# Check Monaco's res genes but correlate with TIMP1
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = monaco.res.set,
         cor.feature = "TIMP1", 
         group.by    = "Patient", 
         ncol        = 5, 
         file.name   = "gene_set_correlations/macrophages.res_genes.TIMP1.pdf",
         width       = 20, 
         height      = 8
)

# Do the id genes
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "CLEC4A", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/macrophages.id_genes.CLEC4A.pdf",
         width       = 20, 
         height      = 20
)

# Do the id genes but correlate them with CD14
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "CD14", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/macrophages.id_genes.CD14.pdf",
         width       = 20, 
         height      = 20
)

# Do the id genes but correlate them with CASP1
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "CASP1", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/macrophages.id_genes.CASP1.pdf",
         width       = 20, 
         height      = 20
)

# Do the id genes but correlate them with TIMP1
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "TIMP1", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/macrophages.id_genes.TIMP1.pdf",
         width       = 20, 
         height      = 20
)

# Do the id genes but correlate them with CD68
plot.cor(object      = final.pop.call.integrated.mye.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "CD68", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/macrophages.id_genes.CD68.pdf",
         width       = 20, 
         height      = 20
)

# Do the id genes but correlate them with KDM5C
plot.cor(object      = final.pop.call.from_full.integrated.mac.seurat, 
         features    = refined.pop.id.genes,
         cor.feature = "KDM5C", 
         group.by    = "Patient", 
         ncol        = 4, 
         file.name   = "gene_set_correlations/macrophages.id_genes.KDM5C.pdf",
         width       = 20, 
         height      = 20
)

## Let's check out the correlations with the mac archetypes
# Set up archetypes
archetypes <- as.vector(Idents(final.pop.call.integrated.mye.seurat))
archetypes[grep("Inflammatory", archetypes)] <- "Inflammatory"
archetypes[grep("Resident", archetypes)]     <- "Resident"
archetypes[grep("Foamy", archetypes)]        <- "Foamy"
archetypes[grep("Lipid", archetypes)]        <- "Resident"

# Sanity check
levels(as.factor(archetypes))

# Add to the seurat object
final.pop.call.integrated.mye.seurat <- AddMetaData(final.pop.call.integrated.mye.seurat,metadata = archetypes, col.name = "archetype")

# Iterate over archetypes
dir.create("gene_set_correlations/archetypes", showWarnings = F, recursive = T)
for(theArchetype in unique(final.pop.call.integrated.mye.seurat$archetype)){
  archetype.tmp.seurat <- subset(final.pop.call.integrated.mye.seurat, subset = archetype == theArchetype)
  cat(c("Working on", theArchetype, "\n"))
  
  # Check Monaco's res genes
  plot.cor(object      = archetype.tmp.seurat, 
           features    = monaco.res.set,
           cor.feature = "CLEC4A", 
           group.by    = "Patient", 
           ncol        = 5, 
           file.name   = paste("gene_set_correlations/archetypes/", theArchetype, ".res_genes.CLEC4A.pdf", sep = ""),
           width       = 20, 
           height      = 8
  )
  
  # Check Monaco's res genes but correlate with CD14
  plot.cor(object      = archetype.tmp.seurat, 
           features    = monaco.res.set,
           cor.feature = "CD14", 
           group.by    = "Patient", 
           ncol        = 5, 
           file.name   = paste("gene_set_correlations/archetypes/", theArchetype, ".res_genes.CD14.pdf", sep = ""),
           width       = 20, 
           height      = 8
  )
  
  # Check Monaco's res genes but correlate with CASP1
  plot.cor(object      = archetype.tmp.seurat, 
           features    = monaco.res.set,
           cor.feature = "CASP1", 
           group.by    = "Patient", 
           ncol        = 5, 
           file.name   = paste("gene_set_correlations/archetypes/", theArchetype, ".res_genes.CASP1.pdf", sep = ""),
           width       = 20, 
           height      = 8
  )
  
  # Check Monaco's res genes but correlate with TIMP1
  plot.cor(object      = archetype.tmp.seurat, 
           features    = monaco.res.set,
           cor.feature = "TIMP1", 
           group.by    = "Patient", 
           ncol        = 5, 
           file.name   = paste("gene_set_correlations/archetypes/", theArchetype, ".res_genes.TIMP1.pdf", sep = ""),
           width       = 20, 
           height      = 8
  )
  
  # Do the id genes
  plot.cor(object      = archetype.tmp.seurat, 
           features    = refined.pop.id.genes,
           cor.feature = "CLEC4A", 
           group.by    = "Patient", 
           ncol        = 4, 
           file.name   = paste("gene_set_correlations/archetypes/", theArchetype, ".id_genes.CLEC4A.pdf", sep = ""),
           width       = 20, 
           height      = 20
  )
  
  # Do the id genes but correlate them with CD14
  plot.cor(object      = archetype.tmp.seurat, 
           features    = refined.pop.id.genes,
           cor.feature = "CD14", 
           group.by    = "Patient", 
           ncol        = 4, 
           file.name   = paste("gene_set_correlations/archetypes/", theArchetype, ".id_genes.CD14.pdf", sep = ""),
           width       = 20, 
           height      = 20
  )
  
  # Do the id genes but correlate them with CASP1
  plot.cor(object      = archetype.tmp.seurat, 
           features    = refined.pop.id.genes,
           cor.feature = "CASP1", 
           group.by    = "Patient", 
           ncol        = 4, 
           file.name   = paste("gene_set_correlations/archetypes/", theArchetype, ".id_genes.CASP1.pdf", sep = ""),
           width       = 20, 
           height      = 20
  )
  
  # Do the id genes but correlate them with TIMP1
  plot.cor(object      = archetype.tmp.seurat, 
           features    = refined.pop.id.genes,
           cor.feature = "TIMP1", 
           group.by    = "Patient", 
           ncol        = 4, 
           file.name   = paste("gene_set_correlations/archetypes/", theArchetype, ".id_genes.TIMP1.pdf", sep = ""),
           width       = 20, 
           height      = 20
  )
}

# Monaco res genes violin plots
customVln(features  = monaco.res.set, 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 5,
          name      = "gene_set_correlations/monaco res markers - violin.pdf", 
          width     = 23, height = 7, draw.names = F, cols = M.int_refined.pop.colors, stack = F)

customVln(features  = "CLEC4A", 
          object    = final.pop.call.integrated.mye.seurat,
          ncol      = 1,
          name      = "gene_set_correlations/CLEC4A - violin.pdf", 
          width     = 7, height = 7, draw.names = F, cols = M.int_refined.pop.colors, stack = F)

