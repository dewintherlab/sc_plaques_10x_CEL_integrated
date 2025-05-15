## Export 10X object for Xenium panel design
final.pop.call.integrated.full.seurat.10X        <- subset(final.pop.call.integrated.full.seurat, Method == "10X")
final.pop.call.integrated.full.seurat.10X.plaque <- subset(final.pop.call.integrated.full.seurat.10X, Tissue == "plaque")

DefaultAssay(final.pop.call.integrated.full.seurat.10X.plaque) <- "RNA"


## Sanity checks
# All integers
all.equal(GetAssayData(object=final.pop.call.integrated.full.seurat.10X.plaque, assay="RNA", slot="counts")@x, as.integer(GetAssayData(object=final.pop.call.integrated.full.seurat.10X.plaque, assay="RNA", slot="counts")@x))
head(GetAssayData(object=final.pop.call.integrated.full.seurat.10X.plaque, assay="RNA", slot="counts")@x)

# Check row names
head(rownames(final.pop.call.integrated.full.seurat.10X.plaque))

# Check meta features
dplyr::glimpse(GetAssay(final.pop.call.integrated.full.seurat.10X.plaque)@meta.features)

# Check metadata
dplyr::glimpse(final.pop.call.integrated.full.seurat.10X.plaque[[]])

# Add cell_types
final.pop.call.integrated.full.seurat.10X.plaque <- AddMetaData(final.pop.call.integrated.full.seurat.10X.plaque, Idents(final.pop.call.integrated.full.seurat.10X.plaque), col.name = "cell_type")

## Add ensembl IDs
symbols <- row.names(GetAssay(final.pop.call.integrated.full.seurat.10X.plaque)@meta.features)
ens     <- AnnotationDbi::select(org.Hs.eg.db, keys = symbols, columns = "ENSEMBL", keytype = "SYMBOL")
d       <- duplicated(ens$SYMBOL)
ens     <- ens[!d,"ENSEMBL"]

tt <- GetAssay(final.pop.call.integrated.full.seurat.10X.plaque)
tt@meta.features$ENSEMBL <- ens

final.pop.call.integrated.full.seurat.10X.plaque[["RNA"]] <- tt

# Load dropletutils to epxort to MEX format
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("DropletUtils", quietly = TRUE))
  BiocManager::install("DropletUtils")

# Load library
library(DropletUtils)

# Run function
write10xCounts(
  "reference_data",
  GetAssayData(final.pop.call.integrated.full.seurat.10X.plaque, assay="RNA", slot="counts"),
  gene.id = GetAssay(final.pop.call.integrated.full.seurat.10X.plaque)@meta.features[["ENSEMBL"]],
  gene.symbol = rownames(final.pop.call.integrated.full.seurat.10X.plaque),
  barcodes = colnames(final.pop.call.integrated.full.seurat.10X.plaque),
  type = "sparse",
  version = "3"
)

# List the files
list.files("reference_data")

# Define function
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {
  
  if (require("data.table", quietly = TRUE)) {
    data.table::fwrite(
      data.table::data.table(
        barcode = barcodes,
        annotation = unlist(data[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv")
    )
  } else {
    write.table(
      data.frame(
        barcode = barcodes,
        annotation = unlist(data[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv"),
      sep = ",", row.names = FALSE
    )
  }
  
  bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))
  
  utils::zip(
    bundle,
    list.files(out_dir, full.names = TRUE),
    zip = "zip"
  )
  
  if (file.info(bundle)$size / 1e6 > 500) {
    warning("The output file is more than 500 MB and will need to be subset further.")
  }
}

# Run function
bundleOutputs(out_dir = "reference_data", data = final.pop.call.integrated.full.seurat.10X.plaque, cell_type = "cell_type")

# Output the seurat object as well
saveRDS(final.pop.call.integrated.full.seurat.10X.plaque, file = "reference_data/scPlaque 10X data.rds")


