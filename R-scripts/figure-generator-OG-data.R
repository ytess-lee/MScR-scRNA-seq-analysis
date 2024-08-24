# figure generator for old RDS

library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SeuratObject)
library(glmGamPoi)

# set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")


integrated_OG <- readRDS("RDS/210524-mechanobiology23_integretatedReps_050424.RDS") 
integrated_OG <- FindNeighbors(integrated_OG, dims = 1:13)
integrated_OG <- FindClusters(integrated_OG, dims= 1:13, resolution = 0.4)

temp_idents <- integrated_OG@active.ident

DimPlot(integrated_OG, group.by = 'seurat_clusters', label = T) +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/OG-data/240804-UMAP-old-dataset-cluster12-478x300.png", width = 478, height = 300)
dev.off()

# assign the culture condition to the active identitites for the DEG analysis
Idents(integrated_OG) <- "Condition" 

# deg <- FindAllMarkers(integrated_OG,
#                        logfc.threshold = 0.25,
#                        min.pct = 0.1,
#                        only.pos = TRUE,
#                        slot = 'counts')

filtered_markers <- deg %>% filter(p_val_adj < 0.05)
# write.csv(filtered_markers, file = 'output-dissertation-figures_240720/OG-data/240713-findallmarker-across-clusters(res0.3)-filtered.csv')

filtered_markers <- read.csv(file = 'output-dissertation-figures_240720/OG-data/240713-findallmarker-across-clusters(res0.3)-filtered-AUTOMATED.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

#write.csv(top10, file = '240713-findallmarker-across-clusters(res0.3)-top10-filtered-temp.csv')


DoHeatmap(integrated_OG, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "output-dissertation-figures_240720/figure2/240804-DEG-Condition-res0.4-top10-omitted-removed-598x767.png", width = 598, height = 767)
dev.off() 

Idents(integrated) <- 'Condition'
DimPlot(integrated_OG, split.by = 'Condition') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/OG-data/240721-umap-seurat-clusters(split-by-Condition)-integrated_OG-478x300.png", width = 478, height = 300)
dev.off()

DimPlot(integrated_OG, split.by = 'Day') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/OG-data/240721-umap-seurat-clusters(split-by-Day)-integrated_OG-478x300.png", width = 478, height = 300)
dev.off()

DimPlot(integrated_OG, split.by = 'dataset') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/OG-data/240804-umap-seurat-clusters(split-by-dataset)-integrated_OG-750x325.png", width = 750, height = 325)
dev.off()
integrated_OG$Condition <- integrated_OG@active.ident

integrated_OG$Condition <- factor(integrated_OG$Condition, levels = c("Prestimulation", "Static", "Pulsatile", "Laminar"))
DimPlot(integrated_OG, group.by = 'Condition') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()+ ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/OG-data/240804-umap-seurat-clusters(by condition)-integrated_OG-478x300.png", width = 478, height = 300)
dev.off()



# subset the eht clusters to check whether we can see that decision point between EHT/ endothelial identity
EHT <- subset(integrated_OG, idents = c('9','8','10','12'))

EHT <- FindClusters(EHT, dims= 1:13, resolution = 0.1)

DimPlot(EHT, group.by = 'seurat_clusters') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/OG-data/240804-UMAP-old-dataset-EHT-297x288.png", width = 297, height = 288)
dev.off()
Idents(EHT)
EHT$seurat_clusters <- factor(EHT$seurat_clusters, levels = c('1','3','2','0','4'))









