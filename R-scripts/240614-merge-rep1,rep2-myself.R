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


####### MERGE the replicates with harmony ########

setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")

# We merge the two datasets.
Rep1<- readRDS("/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/New files without pulsatile harmonised/SeuratObjectMechano_Rep1_for_integration_NOPULSATILE.RDS")
Rep2<- readRDS("/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/New files without pulsatile harmonised/SeuratObjectMechano_Rep2_for_integration_NO_PULSATILE.RDS")

# METHOD 1 - SCTransform before merging
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
saveRDS(Reps, "/Volumes/AF_Seq/Mechanobiology/Integrates_Reps_NO_PULSATILE/240614_mechanobiology23_integretatedReps_NO_PULSATILE_Method1.RDS")

#####################################################################################################################################
# METHOD 2 - NO SCTransform before or after harmony, and try adding SCTransform before merging
# create a list with the normalised objects
norm_seurat_list <- c(Rep1, Rep2)

# Find most variable features across samples to integrate
integ_features <- SelectIntegrationFeatures(object.list = norm_seurat_list, nfeatures = 3000) 

# Merge normalized samples
merged_seurat <- merge(Rep1, Rep2, merge.data = TRUE)

# Set the SCT as the assay for the merged_seurat
DefaultAssay(merged_seurat) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat) <- integ_features

# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)

# We need to specify which cells are from which dataset
# We do this by adding a new column to the meta.data
# We add a column called "dataset", rep repeats a value n times.
# In this case ncol(d1) is the number of columns  (= number of cells) in d1
merged_seurat$dataset <- c(rep("Rep1", ncol(Rep1)),rep("Rep2", ncol(Rep2)))
merged_seurat <- SCTransform(merged_seurat)

# Run Harmony. The group.by.vars argument specifies which column
# in the meta.data to use to group the cells
Reps <- RunHarmony(merged_seurat,
                   group.by.vars = "dataset",
                   reduction = "pca",
                   assay = "SCT",
                   reduction.save = "harmony")

DimPlot(merged_seurat)

Reps <- RunUMAP(Reps, reduction = "harmony", assay = "SCT", dims = 1:40)
Reps <- FindNeighbors(object = Reps, reduction = "harmony")
Reps <- FindClusters(Reps, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2))

saveRDS(Reps, "New files without pulsatile harmonised/merged with Harmony on normalised data_FINAL/240615_pleasework_Mechanobiology_Harmonised_2reps.RDS")

# refer to this tutorial on how to continue
# https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html
# I have followed the second option for the integration

