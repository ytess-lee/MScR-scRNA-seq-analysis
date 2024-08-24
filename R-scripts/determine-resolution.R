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
integrated <- readRDS("RDS/240719-integrated-with-cc-data.RDS") 

integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.3)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis


# merged UAMP
DimPlot(integrated, reduction = "umap")

# res = 0.4 =========================================================================TRUE# res = 0.3 =========================================================================================================================
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.3)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis

DimPlot(integrated, reduction = "umap")
dev.copy(png, file = "240713-CLUSTERS(res0.3)-937x708.png", width = 937, height = 708)
dev.off()  # Close the PNG device

res3 <- FindAllMarkers(integrated,
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             only.pos = TRUE,
                             slot = 'counts')
write.csv(res3, file = '240713-findallmarker-across-clusters(res0.3).csv')
filtered_markers <- res3 %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, file = '240713-findallmarker-across-clusters(res0.3)-filtered.csv')
filtered_markers <- read.csv(file = '/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/output-subset-1-7, 8+9/main(0-5)/240718_filtered(p)_markers_new-data-clusters_endo(0-5)-subset-AUTOMATED.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top10, file = '240713-findallmarker-across-clusters(res0.3)-top10-filtered-temp.csv')

top10 <- read.csv(file = 'output-new/Determining resolution/csv/240804-findallmarker-across-clusters(res0.3)-AUTOMATED.csv')
top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))




DoHeatmap(integrated, features = top10$gene) + theme(axis.text.y = element_text(size = 10))
dev.copy(png, file = "240812-CLUSTERS(res0.3)-top10-temp.png", width = 900, height = 1200)
dev.off()  # Close the PNG device

top20 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(20, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

DoHeatmap(integrated, features = top20$gene) + theme(axis.text.y = element_text(size = 8))


# res = 0.5 =========================================================================================================================

# res = 0.4 =========================================================================TRUE# res = 0.4 =========================================================================================================================
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.4)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis

DimPlot(integrated, reduction = "umap")
dev.copy(png, file = "240702-CLUSTERS(res0.4)-937x708.png", width = 937, height = 708)
dev.off()  # Close the PNG device

deseq.res4 <- FindAllMarkers(integrated,
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             only.pos = TRUE,
                             slot = 'counts')
write.csv(deseq.res4, file = '240702-findallmarker-across-clusters(res0.4).csv')
filtered_markers <- deseq.res4 %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, file = '240702-findallmarker-across-clusters(res0.4)-filtered.csv')
filtered_markers <- read.csv('output-new/Determining resolution/csv/240702-findallmarker-across-clusters(res0.4)-AUTOMATED.csv')
top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

DoHeatmap(integrated, features = top10$gene) + theme(axis.text.y = element_text(size = 10))
dev.copy(png, file = "240702-CLUSTERS(res0.4)-top10-temp-1300x900.png", width = 1300, height = 900)
dev.off()  # Close the PNG device


# res = 0.5 =========================================================================================================================
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.5)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis

DimPlot(integrated, reduction = "umap")
dev.copy(png, file = "240702-CLUSTERS(res0.5)-937x708.png", width = 937, height = 708)
dev.off()  # Close the PNG device

deseq.res4 <- FindAllMarkers(integrated,
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             only.pos = TRUE,
                             slot = 'counts')
write.csv(deseq.res4, file = '240702-findallmarker-across-clusters(res0.5).csv')
filtered_markers <- deseq.res4 %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, file = '240702-findallmarker-across-clusters(res0.5)-filtered.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(top10, file = '240702-findallmarker-across-clusters(res0.5)-top10-filtered-temp.csv')

sum(top10$gene %in% rownames(GetAssayData(integrated, layer = 'scale.data')))
DoHeatmap(integrated, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240702-CLUSTERS(res0.5)-top10-temp-1300x900.png", width = 1300, height = 900)
dev.off()  # Close the PNG device
rm(deseq.res4)
# res = 0.6 =========================================================================================================================
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.6)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis

DimPlot(integrated, reduction = "umap")
dev.copy(png, file = "240702-CLUSTERS(res0.6)-937x708.png", width = 937, height = 708)
dev.off()  # Close the PNG device

deseq.res4 <- FindAllMarkers(integrated,
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             only.pos = TRUE,
                             slot = 'counts')
write.csv(deseq.res4, file = '240702-findallmarker-across-clusters(res0.6).csv')
filtered_markers <- deseq.res4 %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, file = '240702-findallmarker-across-clusters(res0.6)-filtered.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(top10, file = '240702-findallmarker-across-clusters(res0.6)-top10-filtered-temp.csv')

sum(top10$gene %in% rownames(GetAssayData(integrated, layer = 'scale.data')))
DoHeatmap(integrated, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240702-CLUSTERS(res0.6)-top10-temp-1300x900.png", width = 1300, height = 900)
dev.off()  # Close the PNG device

# res = 0.7 =========================================================================================================================
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.7)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis

DimPlot(integrated, reduction = "umap")
dev.copy(png, file = "240702-CLUSTERS(res0.7)-937x708.png", width = 937, height = 708)
dev.off()  # Close the PNG device

deseq.res4 <- FindAllMarkers(integrated,
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             only.pos = TRUE,
                             slot = 'counts')
write.csv(deseq.res4, file = '240702-findallmarker-across-clusters(res0.7).csv')
filtered_markers <- deseq.res4 %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, file = '240702-findallmarker-across-clusters(res0.7)-filtered.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(top10, file = '240702-findallmarker-across-clusters(res0.7)-top10-filtered-temp.csv')

sum(top10$gene %in% rownames(GetAssayData(integrated, layer = 'scale.data')))
DoHeatmap(integrated, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240702-CLUSTERS(res0.7)-top10-temp-1300x900.png", width = 1300, height = 900)
dev.off()  # Close the PNG device