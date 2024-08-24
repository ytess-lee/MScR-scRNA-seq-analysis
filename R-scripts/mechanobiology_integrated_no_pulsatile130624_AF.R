library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(monocle3)
library(DESeq2)
library(harmony)
library(grid)


####### IMPORT DATA TO CREATE SEURAT OBJECT #######

# REPLICATE 1
setwd("/Volumes/AF_Seq/Mechanobiology/Replicate 1 - Jun 2023/CellRanger_output/AllSamples")
# Read the matrix from the CellRanger output
# This is returned as a list with 2 elements:
# `Gene Expression`
# `CRISPR Guide`
data10x <- Read10X(".")

# Import gene expression data into Seurat. This (RNA) will be the default assay 
Rep1 <- CreateSeuratObject(data10x, 
                           min.cells = 3, min.features = 200, 
                           project = "Mechano_Rep1_2023")

# Remove the data that are now loaded into Seurat Object
rm(data10x)

Rep1[["percent.mt"]] <- PercentageFeatureSet(Rep1, pattern = "^MT-")
Rep1 <- subset(Rep1, 
               subset = percent.mt > 0.8 & percent.mt < 9 & 
                 nFeature_RNA > 1700 & nFeature_RNA < 9000)

Rep1[["libraryID"]] <- substr(colnames(Rep1), 18, 18)

# assign the day and treatment ID by creating a new metadata
# ifelse (condition, do if yes, do if no)
Rep1[["Day"]] <- ifelse(Rep1[["libraryID"]] == "1"
                        , "10", "13")
Rep1[["Condition"]] <- "Prestimulation"
Rep1[["Condition"]][Rep1[["libraryID"]] == "2"] <- "Static"
Rep1[["Condition"]][Rep1[["libraryID"]] == "3"] <- "Pulsatile"
Rep1[["Condition"]][Rep1[["libraryID"]] == "4"] <- "Laminar"

# remove pulsatile != means not

Rep1 <- subset(Rep1, Condition != "Pulsatile")

# run the normalisationa and dimension reduction
Rep1 <- SCTransform(Rep1, vars.to.regress = "percent.mt")
Rep1 <- FindVariableFeatures(Rep1)
Rep1 <- RunPCA(Rep1)

# save the RDS file or load the RDS file
saveRDS(Rep1, "/Volumes/AF_Seq/Mechanobiology/Integrates_Reps_NO_PULSATILE/SeuratObjectMechano_Rep1_for_integration_NOPULSATILE.RDS")

# REPLICATE 2
# this follows the same as for Rep1
setwd("/Volumes/AF_Seq/Mechanobiology/Replicate 2 - Jul 2023/Aligned_data/All_samples/outs/count/filtered_feature_bc_matrix")
data10x <- Read10X(".")
Rep2 <- CreateSeuratObject(data10x, 
                           min.cells = 3, min.features = 200, 
                           project = "Mechano_Rep1_2023")
rm(data10x)
Rep2[["percent.mt"]] <- PercentageFeatureSet(Rep2, pattern = "^MT-")
Rep2 <- subset(Rep2, 
               subset = percent.mt > 0.8 & percent.mt < 9 & 
                 nFeature_RNA > 1700 & nFeature_RNA < 9000)

Rep2[["libraryID"]] <- substr(colnames(Rep2), 18, 18)
Rep2[["Day"]] <- ifelse(Rep2[["libraryID"]] == "1"
                        , "10", "13")
Rep2[["Condition"]] <- "Prestimulation"
Rep2[["Condition"]][Rep2[["libraryID"]] == "2"] <- "Static"
Rep2[["Condition"]][Rep2[["libraryID"]] == "3"] <- "Pulsatile"
Rep2[["Condition"]][Rep2[["libraryID"]] == "4"] <- "Laminar"
Rep2 <- SCTransform(Rep2, vars.to.regress = "percent.mt")
Rep2 <- FindVariableFeatures(Rep2)
Rep2 <- RunPCA(Rep2)
Rep2 <- subset(Rep2, Condition != "Pulsatile")
saveRDS(Rep2, "/Volumes/AF_Seq/Mechanobiology/Integrates_Reps_NO_PULSATILE/SeuratObjectMechano_Rep2_for_integration_NO_PULSATILE.RDS")

####### MERGE the replicates with harmony ########

setwd("/Volumes/AF_Seq/Mechanobiology/Integrates_Reps_NO_PULSATILE")
# We merge the two datasets.
Reps <- merge(Rep1, Rep2, merge.data = TRUE)

# We need to specify which cells are from which dataset
# We do this by adding a new column to the meta.data
# We add a column called "dataset", rep repeats a value n times.
# In this case ncol(d1) is the number of columns  (= number of cells) in d1
Reps$dataset <- c(rep("Rep1", ncol(Rep1)),rep("Rep2", ncol(Rep2)))
# Normalize the data
Reps <- SCTransform(Reps)
# remove the files that are now merged
rm(Rep1, Rep2)
# Run PCA
Reps <- RunPCA(Reps)

# Run Harmony. The group.by.vars argument specifies which column
# in the meta.data to use to group the cells
Reps <- RunHarmony(Reps,
                         assay = "SCT",
                         dims.use = 1:30,
                         group.by.vars = "dataset")
saveRDS(Reps, "/Volumes/AF_Seq/Mechanobiology/Integrates_Reps_NO_PULSATILE/mechanobiology23_integretatedReps_NO_PULSATILE130624.RDS")

# = = = = == = = = = = = = = == = = = == = = = = == = = = = = = = = == = = = = = = == = = #
# add those omitted features back into the harmonised dataset AND
# rescale the dataset to include the new features
check <- readRDS('/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/RDS/240615_pleasework_Mechanobiology_Harmonised_2reps.RDS')

# Check variable features
variable_features <- VariableFeatures(check)

# Print or view the variable features to ensure your desired features are included
print(variable_features)

# In DoHeatmap(integrated, features = top20$gene) : The following features were omitted as they were not
# found in the scale.data slot for the SCT assay: 'AL136084.3', 'KLF4', 'SLFN11', 'AL365214.2', 'AC110611.1', 'AL161629.1', 
# 'LAMB3', 'KLHDC8A', 'AC091563.1', 'SERPINB4', 'CAPN8', 'LYPD5', 'CCL8', 'HSPA12B', 'AC016168.4', 'MCEMP1', 'FGR', 'P2RY13', 
# 'SERPINB10', 'BPI', 'LIN28A', 'LINC01629', 'NXPH2', 'AP000692.2', 'AL355075.4', 'VSTM1', 'AL021155.5', 'CCN4', 'Z93241.1', 'AC239799.2'

# Omitted features in downstream analysis
desired_features <- c('AL136084.3', 'KLF4', 'SLFN11', 'AL365214.2', 'AC110611.1', 'AL161629.1',
                      'LAMB3', 'KLHDC8A', 'AC091563.1', 'SERPINB4', 'CAPN8', 'LYPD5', 'CCL8', 'HSPA12B', 'AC016168.4', 
                      'MCEMP1', 'FGR', 'P2RY13', 'SERPINB10', 'BPI', 'LIN28A', 'LINC01629', 
                      'NXPH2', 'AP000692.2', 'AL355075.4', 'VSTM1', 'AL021155.5', 'CCN4', 'Z93241.1', 'AC239799.2' ) # replace with actual gene names

# Check for missing features
missing_features <- setdiff(desired_features, variable_features)

# Add missing features to the variable features if necessary
if (length(missing_features) > 0) {
  VariableFeatures(check) <- unique(c(variable_features, missing_features))
}

Reps <- check

# Rescale the data including the newly added features
Reps <- ScaleData(Reps, features = VariableFeatures(Reps))

# Create the heatmap
DoHeatmap(Reps, features = desired_features) + NoLegend()
saveRDS(Reps, '/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/RDS/240701-add-omitted-no-pulsatile-two-reps-dataset.RDS')

