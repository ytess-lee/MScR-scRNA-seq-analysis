# cluster others: endothelial
# cluster 8: EHT
# cluster 7: blood lineages

# 'Prostate cells','Liver/ Kidney cells','Neural cells','Stromal/ muscle cells','Neural cells','Stromal/ fibroblast cells','Cardiac/ renal cells',
# 'Immune cells','HSC','Immune cells','Neural cells','Neural cells','Erythrocytes'


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
library(Seurat)
library(SeuratObject)
library(glmGamPoi)
library(EnhancedVolcano)
library(patchwork)

# set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")

# using new RDS file from Anto - 240614
integrated <- readRDS("RDS/240615_pleasework_Mechanobiology_Harmonised_2reps.RDS") 

integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.5)
DimPlot(integrated, reduction = "umap", label = TRUE, label.size = 5)


# ======================== # ANNOTATION - make res0.5 to this # ======================== # 
# Define new cluster IDs
new.cluster.ids <- c('0-Prostate cells','1-Liver/ kidney cells','2-Neural cells','3-Stromal/ muscle cells',
                     '4-Neural cells','5-Stromal/ fibroblast cells','6-Cardiac cells','7-Haematopoietic lineages',
                     '8-Haematopoietic stem cells','9-Immune cells','10-Neural cells','11-Neural cells',
                     '12-Erythrocytes')

# Ensure the new cluster IDs have names that match current cluster levels
current_clusters <- levels(integrated@meta.data$SCT_snn_res.0.5)
if (length(current_clusters) != length(new.cluster.ids)) {
  stop("Number of new cluster IDs does not match the number of current clusters")
}
names(new.cluster.ids) <- current_clusters

# Map the new cluster IDs to the existing SCT_snn_res.0.5 column
integrated@meta.data$SCT_snn_res.0.5 <- new.cluster.ids[as.character(integrated@meta.data$SCT_snn_res.0.5)]

# Convert to factor with proper levels
integrated@meta.data$SCT_snn_res.0.5 <- factor(integrated@meta.data$SCT_snn_res.0.5, levels = unique(new.cluster.ids))

clusters_12 <- DimPlot(integrated, reduction = 'umap', group.by = "SCT_snn_res.0.5", label = TRUE, pt.size = 0.5, label.size=5) +NoLegend()
clusters_12

# ======================== # ANNOTATION make seurat_clusters# ======================== # 
# Define new cluster IDs
new.cluster.ids.endo <- c('Endothelial cells','Endothelial cells','Endothelial cells','Endothelial cells',
                     'Endothelial cells','Endothelial cells','Endothelial cells','7-Haematopoietic lineages',
                     '8-Haematopoietic stem cells','9-Immune cells','Endothelial cells','Endothelial cells',
                     '12-Erythrocytes')

# Ensure the new cluster IDs have names that match current cluster levels
current_clusters_endo <- levels(integrated@meta.data$seurat_clusters)
if (length(current_clusters_endo) != length(new.cluster.ids.endo)) {
  stop("Number of new cluster IDs does not match the number of current clusters")
}
names(new.cluster.ids.endo) <- current_clusters_endo

# Map the new cluster IDs to the existing SCT_snn_res.0.5 column
integrated@meta.data$seurat_clusters <- new.cluster.ids.endo[as.character(integrated@meta.data$seurat_clusters)]

# Convert to factor with proper levels
integrated@meta.data$seurat_clusters <- factor(integrated@meta.data$seurat_clusters, levels = unique(new.cluster.ids.endo))

by_GO <- DimPlot(integrated, reduction = "umap",group.by = 'seurat_clusters', label = TRUE, pt.size = 0.5, label.size=5) + NoLegend()
by_GO

clusters_12|by_GO
dev.copy(png, file = "240707-side-by-side-annotation-GO+cell-marker-UMAP-1059x708-new.png", width = 1059, height = 800)
dev.off()


# for subset 6+7
# 0	Endothelial cell priming for EHT	
# 1	Endothelial cell	
# 2	EHT	
# 3	Myleiod linage + multi-potential progenitors	
# 4	Lymphoid lineage	


# ======================== # ANNOTATION - make res0.5 to this # ======================== # 
# Define new cluster IDs
cluster6_7 <- readRDS(file = 'RDS/240714-cluster6+7(EHT+Blood)-subset.RDS')
cluster6_7 <- FindClusters(cluster6_7, dims= 1:15, resolution = 0.2)
DimPlot(cluster6_7, group.by = 'SCT_snn_res.0.2', label = T)
new.cluster.ids <- c('Endothelial cells priming for EHT', 
                     'Myeloid lineage + MPP',
                     'Lymphoid lineage',
                     'EHT', 
                     'Endothelial cells')

# Ensure the new cluster IDs have names that match current cluster levels
current_clusters <- levels(cluster6_7@meta.data$SCT_snn_res.0.2)
if (length(current_clusters) != length(new.cluster.ids)) {
  stop("Number of new cluster IDs does not match the number of current clusters")
}
names(new.cluster.ids) <- current_clusters

# Map the new cluster IDs to the existing SCT_snn_res.0.5 column
cluster6_7@meta.data$SCT_snn_res.0.2 <- new.cluster.ids[as.character(cluster6_7@meta.data$SCT_snn_res.0.2)]

# Convert to factor with proper levels
#cluster6_7@meta.data$SCT_snn_res.0.2 <- factor(cluster6_7@meta.data$SCT_snn_res.0.5, levels = unique(new.cluster.ids))

clusters_67 <- DimPlot(cluster6_7, reduction = 'umap', group.by = "SCT_snn_res.0.2", label = TRUE, pt.size = 0.5, label.size=5) +NoLegend()
clusters_67
DimPlot(cluster6_7, group.by = 'seurat_clusters', label = T)
# ======================== # ANNOTATION make seurat_clusters# ======================== # 
# Define new cluster IDs
DimPlot(integrated, group.by = 'SCT_snn_res.0.3', label = T)
new.cluster.ids.endo <- c('0-Endothelial cells',
                          '1-Endothelial cells',
                          '2-Endothelial cells',
                          '3-Endothelial cells',
                          '4-Endothelial cells',
                          '5-Endothelial cells',
                          '6-Haematopoietic lineage',
                          '7-EHT',
                          '8-Megakaryocytes',
                          '9-Erythroid')

# Ensure the new cluster IDs have names that match current cluster levels
current_clusters_endo <- levels(integrated@meta.data$seurat_clusters)
if (length(current_clusters_endo) != length(new.cluster.ids.endo)) {
  stop("Number of new cluster IDs does not match the number of current clusters")
}
names(new.cluster.ids.endo) <- current_clusters_endo

# Map the new cluster IDs to the existing SCT_snn_res.0.5 column
integrated@meta.data$seurat_clusters <- new.cluster.ids.endo[as.character(integrated@meta.data$seurat_clusters)]

# Convert to factor with proper levels
integrated@meta.data$seurat_clusters <- factor(integrated@meta.data$seurat_clusters, levels = unique(new.cluster.ids.endo))

by_GO <- DimPlot(integrated, reduction = "umap",group.by = 'seurat_clusters', label = TRUE, pt.size = 0.5, label.size=5) + NoLegend()
by_GO

clusters_67|by_GO
dev.copy(png, file = "2407014-cluster6+7-integrated-annotation-GO+cell-marker-UMAP-1500x708-new.png", width = 1500, height = 800)
dev.off()


DimPlot(integrated, group.by = 'SCT_snn_res.0.3', label = 5, pt.size = 0.5, label.size = 6)
dev.copy(png, file = '240714-integrated-9clusters-res0.3-1059x800.png', width = 1059, height = 800)
dev.off()
