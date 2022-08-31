##==========================================================================================================
## Overlay scATAC data
##=========================================================================================================
dir.create("scATAC_integration", showWarnings = F)

## Load data
full.scATAC.signac <- readRDS(file = "raw_data/Kaikkonen_lab_athero_scatac.SignacV1.harmony_2021_12_10.RDS")
full.scATAC.signac


## Use our scRNA-seq data as reference
# First transform the ATAC data with SCT to match the RNN reference
full.scATAC.signac <- SCTransform(full.scATAC.signac, assay = "peaksMACS2")
