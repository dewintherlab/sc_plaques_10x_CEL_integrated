#=================================================================================================================================
## Overview of endpoint or branchpoint objects, saved to disk.
## Can be used to restart analysis at a midway or end point without having to re-run the whole thing.
##================================================================================================================================
dir.create("Seurat_Objects", showWarnings = F)

#=================================================================================================================================
## Objects to save (commented out so we don't accidentally overwrite our saved stuff)
##================================================================================================================================
## Gene ontology data
# saveRDS(human, "Seurat_Objects/human.mart.data.RDS")
# saveRDS(mouse, "Seurat_Objects/mouse.mart.data.RDS")
# 
# ## Seurat objects
# # List of 6 raw seurat objects, freshly loaded
# saveRDS(samples.seurat, "Seurat_Objects/main.seurat.raw.RDS")
# 
# # List of 3 seurat objects (per patient), SCtransformed, embedded, % mt etc. calculated
# saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.prefilter.RDS")
# 
# # List of 3 seurat objects (per patients), SCtransformed, embedded, % mt etc. filtered
# saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.filtered.RDS")
# 
# # Integrated seurat object, embedded
# saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.filtered.integrated.RDS")
# 
# # Integrated seurat object, embedded, mapped to reference
# saveRDS(patients.seurat, "Seurat_Objects/main.patients.seurat.filtered.integrated.mapped.RDS")
#
# # Integrated seurat object, embedded, mapped to reference, populations refined
# saveRDS(refined.patients.seurat, "Seurat_Objects/main.patient.integrated.mapped.refined_populations.seurat.RDS")
#
# # Myeloid cells seurat object
# saveRDS(mye.patients.seurat, "Seurat_Objects/myeloid.cells.seurat.RDS")
# 
# # PBMC only monocytes eurat object
# saveRDS(pbmc.monocytes.seurat, "Seurat_Objects/pbmc.monocytes.seurat.RDS")
#
# # Plaque only macrophages seurat object
# saveRDS(plaque.macrophages.seurat, "Seurat_Objects/plaque.macrophages.seurat.RDS")
#
# ## Marker genes
# # All markers from all patients and idents
# saveRDS(all.patients.seurat.markers, "Seurat_Objects/main.patients.markers.RDS")
# 
# # PBMC monocyte markers
# saveRDS(pbmc.monocytes.seurat.markers, "Seurat_Objects/PBMC.monocytes.markers.RDS")
#
# # Plaque macrophage markers
# saveRDS(plaque.macrophages.seurat.markers, "Seurat_Objects/plaque.macrophages.markers.RDS")
# 
# # Myeloid cells refined Seurat object
# saveRDS(mye.patients.seurat, "Seurat_Objects/myeloid.cells.refined.seurat.RDS")
#
# # Integrated 43 patietns cel-seq and 10X myeloid cells seurat object
# saveRDS(integrated.mye.seurat, "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")
#
# # Full set 43 patients cel-seq and 10X cells seurat object
# saveRDS(integrated.full.seurat, "Seurat_Objects/full.43p_10X.integrated.seurat.RDS")
#
# # ssGSEA data ran on mye.integrated.seurat
# saveRDS(ssGSEA.integrated.mye.seurat, "Seurat_Objects/ssGSEA_data.43p_10X.integrated.RDS")
#
#
# ## Gene Ontologies
# 
# 
# ## Monocle objects
# 
# 
# ## Control samples
# saveRDS(ctrl.combined,         "Seurat_Objects/control.combined.Wang.seurat.RDS")
# saveRDS(ctrl_patient.combined, "Seurat_Objects/control_patient.combined.Wang.seurat.RDS")
# 
# ## Random odds and ends
# saveRDS(vars.to.regress, "Seurat_Objects/vars.to.regress.Rds")


#=================================================================================================================================
## Objects to load
##================================================================================================================================
## Gene ontology data
ont_sets <- readList("cp_and_h.gmt.txt")
human    <- readRDS(file = "Seurat_Objects/human.mart.data.RDS")
mouse    <- readRDS(file = "Seurat_Objects/mouse.mart.data.RDS")

## Seurat objects
# List of 6 raw seurat objects, freshly loaded
samples.seurat <- readRDS(file = "Seurat_Objects/main.seurat.raw.RDS")

# List of 3 seurat objects (per patient), SCtransformed, embedded, % mt etc. calculated
patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.prefilter.RDS")

# List of 3 seurat objects (per patient), SCtransformed, embedded, % mt etc. filtered
patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.filtered.RDS")

# Integrated seurat object, embedded
patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.integrated.RDS")

# Integrated seurat object, embedded, mapped to reference
patients.seurat <- readRDS(file = "Seurat_Objects/main.patients.seurat.filtered.mapped.RDS")

# Integrated seurat object, embedded, mapped to reference, populations refined
refined.patients.seurat <- readRDS(file = "Seurat_Objects/main.patient.integrated.mapped.refined_populations.seurat.RDS")

# Myeloid cells seurat object
mye.patients.seurat <- readRDS(file = "Seurat_Objects/myeloid.cells.seurat.RDS")

# PBMC only monocytes seurat object
pbmc.monocytes.seurat <- readRDS(file = "Seurat_Objects/pbmc.monocytes.seurat.RDS")

# Plaque only macrophages seurat object
plaque.macrophages.seurat <- readRDS(file = "Seurat_Objects/plaque.macrophages.seurat.RDS")

# Myeloid cells refined Seurat object
mye.patients.seurat <- readRDS(file = "Seurat_Objects/myeloid.cells.refined.seurat.RDS")

# Integrated 43 patietns cel-seq and 10X myeloid cells seurat object
integrated.mye.seurat <- readRDS(file = "Seurat_Objects/mye.43p_10X.integrated.seurat.RDS")

# Full set 43 patients cel-seq and 10X cells seurat object
integrated.full.seurat <- readRDS(file = "Seurat_Objects/full.43p_10X.integrated.seurat.RDS")

# ssGSEA data ran on mye.integrated.seurat
ssGSEA.integrated.mye.seurat <- readRDS(file = "Seurat_Objects/ssGSEA_data.43p_10X.integrated.RDS")


## Marker genes
# All markers from all patients and idents
all.patients.seurat.markers    <- readRDS(file = "Seurat_Objects/main.patients.markers.RDS")
pbmc.monocytes.seurat.markers <- readRDS(file = "Seurat_Objects/PBMC.monocytes.markers.RDS")
plaque.macrophages.seurat.markers <- readRDS(file = "Seurat_Objects/plaque.macrophages.markers.RDS")

## Gene Ontologies


## Monocle objects


## Control samples
ctrl.combined         <- readRDS(file = "Seurat_Objects/control.combined.Wang.seurat.RDS")
ctrl_patient.combined <- readRDS(file = "Seurat_Objects/control_patient.combined.Wang.seurat.RDS")

## Random odds and ends
vars.to.regress <- readRDS(file = "Seurat_Objects/vars.to.regress.Rds")