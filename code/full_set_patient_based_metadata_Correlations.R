#=========================================================================================================================
## Analyze correlations with metadata
##========================================================================================================================
dir.create("full_set__patient_based_metadata_correlations/Individual trait distribution", showWarnings = F, recursive = T)
dir.create("full_set__patient_based_metadata_correlations/Individual trait distribution binned", showWarnings = F, recursive = T)


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
# Fetch patient information per cell
md.df <- data.frame(Patient=full.43p.seurat$Patient)

# Expand metadata to all cells
md.df <- merge(md.df, meta.data, by.x = 1, by.y = 2, all.x = T)
md.df$Patient <- NULL

# And add to the Seurat object
for (i in names(md.df)){
  full.43p.seurat <- AddMetaData(full.43p.seurat, md.df[,i], col.name = i)
}

##=================================================================================================================================
## Get per patient distribution of traits

patient.trait.dist <- list()
for(thePatient in unique(full.43p.seurat$Patient)){
  cat(paste("Working on patient", thePatient,"\n\n"))
  # Extract patient cells
  tmp.patient.seurat <- subset(full.43p.seurat, Patient == thePatient)
  
  ## Plot individual relative distributions per trait
  # Get distributions per trait, per population
  trait.dist <- list()
  i          <- 1
  i.max      <- length(colnames(meta.data[-2])) #-2 because we merged on that column
  for(theTrait in colnames(meta.data[-2])){
    cat(paste("Working on trait: ", theTrait, " (", i, " of ", i.max, ")\n", sep =""))
    
    # make a df in the list and populate with the distribution for this trait over all cells
    trait.dist[[theTrait]]            <- data.frame(table(tmp.patient.seurat@meta.data[,theTrait]))
    row.names(trait.dist[[theTrait]]) <-  trait.dist[[theTrait]]$Var1
    trait.dist[[theTrait]]$Var1       <- NULL
    colnames(trait.dist[[theTrait]])  <- "Total"
    
    # Then add the distribution per population
    for(thePop in unique(Idents(tmp.patient.seurat))){
      cat(paste("\tAdding distribution for pop: ", thePop, "...\n", sep =""))
      
      # Subset out the desired population and get the distribution
      tmp.pop.seurat        <- subset(tmp.patient.seurat, idents = thePop)
      tmp.pop.df            <- data.frame(table(tmp.pop.seurat@meta.data[,theTrait]))
      row.names(tmp.pop.df) <- tmp.pop.df$Var1
      tmp.pop.df$Var1       <- NULL
      colnames(tmp.pop.df)  <- thePop
      
      ## Add this pop's distribution to the overall df in the list
      # We need to merge instead of cbind just in case a pop doesn't have a certain category of the trait.
      trait.dist[[theTrait]]            <- merge(trait.dist[[theTrait]], tmp.pop.df, by.x = 0, by.y = 0, all.x = T)
      row.names(trait.dist[[theTrait]]) <- trait.dist[[theTrait]]$Row.names
      trait.dist[[theTrait]]$Row.names  <- NULL
    }
    
    # Get rid of NAs introduced by trait cats not in certain populations
    trait.dist[[theTrait]][is.na(trait.dist[[theTrait]])] <- 0
    
    i <- i + 1
  }
  head(trait.dist)
  
  # Clean up: remove traits with one observation (shouldn't be there to begin with)
  for(theTrait in names(trait.dist)){
    if(dim(trait.dist[[theTrait]])[1] < 2){
      cat(paste("Removed trait: ", theTrait, "\n", sep = ""))
      trait.dist[[theTrait]] <- NULL
    }
  }
  
  # Add the trait dist to the main patient.trait.dist list
  patient.trait.dist[[thePatient]] <- trait.dist
}

names(patient.trait.dist)
patient.trait.dist[["4440"]]


## Add population info to the metadata
# Link populations to idents (uniquely)
ident.table           <- merge(full.43p.seurat@meta.data[,c("Patient", "celltypes.43p"), drop = F],  as.data.frame(Idents(full.43p.seurat)), by = 0)
d                     <- duplicated(ident.table$Patient)
ident.table           <- ident.table[!d,]
ident.table$Row.names <- NULL
colnames(ident.table) <- c("Patient", "Celltype", "Population")

# Merge it with the meta data
pop.meta.data <- merge(meta.data, ident.table, by.x = 2, by.y = 1)
head(pop.meta.data)

pop.meta.data[]   <- lapply(pop.meta.data, factor)
num.pop.meta.data <- sapply(pop.meta.data, unclass)
num.pop.meta.data
chisq.test(num.pop.meta.data)


