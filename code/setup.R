## Install packages
# install.packages(c("BiocManager", "devtools", "rmarkdown", "Rcpp", "RcppArmadillo", "xlsx))
# BiocManager::install(c("Seurat", "ggplot2", "patchwork", "dplyr", "reshape2", "pheatmap", "varhandle", "ggsci", "limma", "ggpubr", "AnnotationDbi", "org.Hs.eg.db", "msigdbr", "dittoSeq", "EGSEA", "biomaRt", "escape", "IRanges", "gage", "pcaMethods", "ggsignif","qusage"))
# devtools::install_github("mojaveazure/seurat-disk")
# devtools::install_github("velocyto-team/velocyto.R")
# remotes::install_github('satijalab/seurat-wrappers')

## Monocle 3:
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')
#
# renv::snapshot()


## Load Packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ggsci)
library(SeuratDisk)
library(limma)
library(ggpubr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EGSEA)
library(biomaRt)
library(velocyto.R)
library(monocle3)
library(SeuratWrappers)
library(ggsignif)
library(xlsx)
library(qusage)
library(escape)
library(msigdbr)
library(dittoSeq)
library(rstatix)
library(varhandle)

# Load custom functions
source("code/functions.R")

#Set up gene ontology objects for pathway analyses
ont_sets <- readList("cp_and_h.gmt.txt")

#Initialize bioMarts
marts <- listMarts()
mart <- marts[1,1]

human <- useMart(mart, dataset = "hsapiens_gene_ensembl")
mouse <- useMart(mart, dataset = "mmusculus_gene_ensembl")

## Order of execution
# setup.R
# load_data.R
# normalise_and_filter.R
# integrate.R
# auto_map_cell_types.R
# DE_Calling_and_cluster_annotation.R
# Myeloid_cells.R
# PBMC_monocytes.R
# plaque_macrophages.R
# myeloid_Cell_identity_back_mapping.R