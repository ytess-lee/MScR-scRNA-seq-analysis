# This script is to identify omitted features, remove them from the exported csv file, reload the file and create top10/20 heatmap
# cell typing based on marker distribution

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

# temporarly store the cluster idents to be able to re-use them later
temp_idents <- integrated@active.ident

# assign the culture condition to the active identitites for the DEG analysis
Idents(integrated) <- "Condition" 

# Reorder the levels of the condition factor
integrated$Condition <- factor(integrated$Condition, levels = c("Prestimulation", "Static", "Laminar"))

# determine the # of dimensions/ PCs we need to capture the majority of the variation in the integrated
ElbowPlot(integrated, ndims = 30)

# # # # # 
#  UMAP #
# # # # # 
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.5)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis


# merged UAMP
DimPlot(integrated, reduction = "umap")
# split UMAP by conditions
DimPlot(integrated, reduction = "umap",
        split.by = "dataset", group.by = "Condition",ncol = 2) # number of columns for display when combining plots


# # # # # # # # # # # # # # # # # # # # ## # # # # 
# FILTERED MARKERS W/ p_val_adj < 0.05 && EXPORT #
# # # # # # # # # # # # # # # # # # # # ## # # # # 
# Apply SCTransform normalization
# integrated <- SCTransform(integrated, verbose = TRUE)
# 
# # Prepare SCT object for marker detection
# integrated <- PrepSCTFindMarkers(integrated)

filtered_markers <- read.csv(file = '240614_filtered(p)_markers-omitted-features-removed-0704.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) 

write.csv(top10, "240614_filtered(p)_markers_top10-omitted-features-removed-0704.csv", row.names = FALSE)
DoHeatmap(integrated, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240704-top10-omitted-features-removed-0704-1059x705.png", width = 1059, height = 705)
dev.off()  # Close the PNG device

top20 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(20, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 
write.csv(top20, "240614_filtered(p)_markers_top20-omitted-features-removed-0704.csv", row.names = FALSE)
DoHeatmap(integrated, features = top20$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240704-top20-omitted-features-removed-0704-1059x705.png", width = 1059, height = 705)
dev.off()  # Close the PNG device


# encounter issue with SCT assay, try mapping the markers from AF's csv
# interestingy it still omitted quite a lot markers i guess it's a no go
AF_marker <- read.csv(file='Marker_Clusters_Mechano_merged.csv')
top10 <- AF_marker %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 
DoHeatmap(integrated, features = top10$gene) + theme(axis.text.y = element_text(size = 8))


# check feature expression
# Prestim: 'GAL','CKS1B','AC012447.1','TEX14','ADAMTS18','SOST','DMKN','ADAMTS20','REN','CD3D'
FeaturePlot(object = integrated, features = c('GAL','CKS1B','AC012447.1','TEX14','ADAMTS18','SOST','DMKN','ADAMTS20','REN','CD3D'))
dev.copy(png, file = "240704-top10-expression-prestim(omitted-removed)-900x700.png", width = 900, height = 700)
dev.off() 

DotPlot(object = integrated, features = c('GAL','CKS1B','AC012447.1','TEX14','ADAMTS18','SOST','DMKN','ADAMTS20','REN','CD3D'), group.by = "Condition") +
  ggtitle("TOP10 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240630-top10-expression-prestim(omitted-removed)-dotplot-1059x705.png", width = 1059, height = 705)
dev.off() 

# Static: 'HLA-DRA','RETN','MS4A4E','IL1RN','CD1E','MS4A7','CD1C','CD14','HLA-DQA1','ELANE','DEFA4'
FeaturePlot(object = integrated, features = c('HLA-DRA','RETN','MS4A4E','IL1RN','CD1E','MS4A7','CD1C','CD14','HLA-DQA1','ELANE','DEFA4'))
dev.copy(png, file = "240704-top10-expression-static(omitted-removed)-900x700.png", width = 900, height = 700)
dev.off() 

DotPlot(object = integrated, features = c('HLA-DRA','RETN','MS4A4E','IL1RN','CD1E','MS4A7','CD1C','CD14','HLA-DQA1','ELANE','DEFA4'), group.by = "Condition") +
  ggtitle("TOP10 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240630-top10-expression-static(omitted-removed)-dotplot-1059x705.png", width = 1059, height = 705)
dev.off() 

# Laminar:'KLF2','LAMA4','COL15A1','NCALD','A2M','FAM167B','UNC13A','IGFBP5','FABP3','DEPTOR'
FeaturePlot(object = integrated, features = c('KLF2','LAMA4','COL15A1','NCALD','A2M','FAM167B','UNC13A','IGFBP5','FABP3','DEPTOR'))
dev.copy(png, file = "240704-top10-expression-laminar(omitted-removed)-900x700.png", width = 900, height = 700)
dev.off() 

DotPlot(object = integrated, features = c('KLF2','LAMA4','COL15A1','NCALD','A2M','FAM167B','UNC13A','IGFBP5','FABP3','DEPTOR'), group.by = "Condition") +
  ggtitle("TOP10 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240630-top10-expression-laminar(omitted-removed)-dotplot-1059x705.png", width = 1059, height = 705)
dev.off() 

# # # BY CELL TYPE

# HSC: 'RUNX1','MLLT3', 'HOXA9', 'MECOM', 'HLF', 'SPINK2'
# CD93: strongly expressed in platelets and megakaryocytes, endothelial, NK and monocytes
# ADGRG1: enriched for functional hHSCs
# MADCAM1:  endothelial cell adhesion molecule - too disperse -> removed
# Not found: CD41(ITGA2B!!), CD43, CD45, kit, CD51, CD61, SCA1, SCA-1 , GPR56, ACRL1, LKZF2, LKZF1, GM28112, HLF,  AL467606 
FeaturePlot(integrated, features = c('CD44', 'RUNX1', 'CD93','GFI1', 'ADGRG1'))
dev.copy(png, file = "240704-HSC-MARKERS-900X700.png", width = 900, height = 700)
dev.off() 

# ERY: 'MLLT3', 'GMPR'
# MLLT3: a regulator of early erythroid and megakaryocytic cell fate
# not found: C2PRF88
FeaturePlot(integrated, features = c('MLLT3','GMPR','KLF1','ITGA2B'))
dev.copy(png, file = "240704-Erythrocytes-MARKERS-900X700.png", width = 900, height = 700)
dev.off() 

# EHT: CD44

# HSC + ENDO: 'CD34', 'PECAM1', 'CDH5', 'CALCRL', 'ROCR', 'ESAM', 'TIE1', 'EMCN'
FeaturePlot(integrated, features = c('CD34', 'PECAM1', 'CDH5', 'CALCRL', 'ESAM', 'TIE1', 'EMCN'))
dev.copy(png, file = "240704-HSC+Endo-MARKERS-900X700.png", width = 900, height = 700)
dev.off() 

# ENDOTHELIAL MARKERS
# NOT FOUND: GPB4
FeaturePlot(integrated, features = c('TEK','PCDH12','PTPRB' , 'RAMP2', 'PECAM1', 'KDR',  'CDH5', 'COL4A2' ,'SOX7','KRT18' ,'EMB' ))
dev.copy(png, file = "240704-ENDOTHELIAL-MARKERS-1-900X700.png", width = 900, height = 700)
dev.off() 

# ENDO:  'RAB38'
FeaturePlot(integrated, features = c('NFE2','DOK2','LGALS9', 'HOPX','MECOM', 'PTK2B', 'MCTP1','MYB', 'CORO1A','SPI1'))            
dev.copy(png, file = "240704-Endothelial-MARKERS-2-900X700.png", width = 900, height = 700)
dev.off()     

# ======================== # BY LOCATION DISTRIBUTION # ======================== # 
# BOTTLE NECK
# remove GMPR: its a GMP reductase and I couldn't find literature supporting its related to EHT
FeaturePlot(integrated, features = c('GATA2', 'ITGA2B', 'CD44' ,'RUNX1' ,'CD82','DOK2', 'LGALS9', 'PTK2B','MCTP1', 'CORO1A', 'SPI1' ))
dev.copy(png, file = "240704-bottleneck-MARKERS-900X700.png", width = 900, height = 700)
dev.off() 

# BLOOD LINEAGES
# ery: 'KLF1', 'GATA1', 'NFE2', 'GFI1'
# myeloid: 'SPI1', 'MYB', 'GATA2'
# remove MCTP1 and LGALS9: couldn't find a direct link with blood lineages
# GATA1 looks funny
FeaturePlot(integrated, features = c('NFE2', 'DOK2', 'PTK2B', 'CORO1A', 'SPI1' , 'MYB', 'KLF1', 'GATA1', 'GFI1','TRPV4' ))            
dev.copy(png, file = "240704-Blood-lineags-900X700.png", width = 900, height = 700)
dev.off()     


# ENDOTHELIAL(B4 BOTTLE NECK)
FeaturePlot(integrated, features = c('TEK','PCDH12', 'RAMP2', 'PECAM1','HOPX', 'MECOM','KDR',  'CDH5', 'COL4A2','KRT18'))
dev.copy(png, file = "240704-Endothelial-B4-EHT-900X700.png", width = 900, height = 700)
dev.off()

# Calcium binding activity
FeaturePlot(integrated, features = c('MCTP1', 'MCTP2', 'TRPV4'))



# dotplot to visualise the distribution of gene experssion across clusters
DotPlot(object = integrated, features = c('PIEZO1', 'PIEZO2', 'AQP1','GATA2', 'ITGA2B', 'CD44' ,'RUNX1' ,'CD82','NFE2', 'DOK2', 'PTK2B', 'CORO1A', 'SPI1' , 'MYB', 'KLF1', 'GATA1', 'GFI1' ), group.by = "Condition") +
  ggtitle("TOP10 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240706-Transporters+EHT+Bottleneck-dotplot-900x705.png", width = 900, height = 705)
dev.off() 

