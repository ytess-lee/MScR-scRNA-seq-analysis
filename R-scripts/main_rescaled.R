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

# set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")

# using new RDS file from Anto - 240614
integrated <- readRDS("RDS/240615_pleasework_Mechanobiology_Harmonised_2reps.RDS") 

# using re-scaled data 
recaled <- readRDS("RDS/240701-add-omitted-no-pulsatile-two-reps-dataset.RDS") 
integrated <- recaled

# # # # # # # # # # # # # # # # #
# DE GENES GROUPED BY CLUSTERS #
# # # # # # # # # # # # # # # # #
# temporarly store the cluster idents to be able to re-use them later
# object@slotName: to access a specific slot within an S4 object (OOP, not seurat specific)
temp_idents <- integrated@active.ident

# assign the culture condition to the active identitites for the DEG analysis
Idents(integrated) <- "Condition" 

# Reorder the levels of the condition factor
integrated$Condition <- factor(integrated$Condition, levels = c("Prestimulation", "Static", "Laminar"))

# determine the # of dimensions/ PCs we need to capture the majority of the variation in the integrated
ElbowPlot(integrated, ndims = 20)


# # # # # 
#  UMAP #
# # # # # 
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.5)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis


# merged UAMP
DimPlot(integrated, reduction = "umap")
dev.copy(png, file = "240701-rescaled-UMAP-merged-973x708.png", width = 973, height = 708)
dev.off()  # Close the PNG device
# split UMAP by conditions
DimPlot(integrated, reduction = "umap",
        split.by = "dataset", group.by = "Condition",ncol = 2) # number of columns for display when combining plots
dev.copy(png, file = "240701-rescaled-UMAP-973x708.png", width = 973, height = 708)
dev.off()  # Close the PNG device

# # # # # 
# tSNE  # - 240606 Didn't save tSNE coz it doesn't say much about the clusters and it's ugly lol
# # # # #
# data <- FindNeighbors(object = integrated, dims = 1:13, verbose = FALSE)
tSNE <- RunTSNE(integrated, dims = 1:13)
# merged tSNE
DimPlot(tSNE, reduction ="tsne")
dev.copy(png, file = "240701-rescaled-tSNE-merged-973x708.png", width = 973, height = 708)
dev.off()  # Close the PNG device
# split tSNE by conditions
DimPlot(tSNE, reduction = "tsne",
        split.by = "dataset", group.by = "Condition", ncol = 2)
dev.copy(png, file = "240701-rescaled-tSNE-973x708.png", width = 973, height = 708)
dev.off()  # Close the PNG device
# # # # # #
# TABLES  #
# # # # # #
table(integrated@active.ident)
table(Idents(integrated))


# # # # # # # # # # # # # # # # # # # # ## # # # # # # # # #
# # Check gene marker expression distribution in clusters #
# # # # # # # # # # # # # # # # # # # # ## # # # # # # # # 
# # Prestimulation: 'SOST', 'REN', 'AC239799.2', 'Z93241.1', 'AC012447.1', 'CCN4', 'GAL', 'ADAMTS18', 'AL021155.5', 'VSTM1'
# Top10, only GAL express something
FeaturePlot(object = integrated, features = c('SOST', 'REN', 'AC239799.2', 'Z93241.1', 'AC012447.1', 'CCN4', 'GAL', 'ADAMTS18', 'AL021155.5', 'VSTM1'))
dev.copy(png, file = "240701-top10-expression-prestim-1059x705.png", width = 1059, height = 705)
dev.off()  # Close the PNG device

DotPlot(object = integrated, features = c('SOST', 'REN', 'AC239799.2', 'Z93241.1', 'AC012447.1', 'CCN4', 'GAL', 'ADAMTS18', 'AL021155.5', 'VSTM1'), group.by = "Condition") +
  ggtitle("TOP10 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240701-top10-expression-prestim-dotplot-900x700.png", width = 900, height = 700)
dev.off()  # Close the PNG device

# Top20: High CKS1B but leaked to the static cluster
FeaturePlot(object = integrated, features = c('AL355075.4','AP000692.2','TEX14','ADAMTS20','NXPH2','CD3D','LINC01629','CKS1B','LIN28A','DMKN'))
dev.copy(png, file = "240701-top20-expression-prestim-1059x705.png", width = 1059, height = 705)
dev.off()  # Close the PNG device

DotPlot(object = integrated, features = c('AL355075.4','AP000692.2','TEX14','ADAMTS20','NXPH2','CD3D','LINC01629','CKS1B','LIN28A','DMKN'), group.by = "Condition") +
  ggtitle("TOP20 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240701-top20-expression-prestim-dotplot-900x700.png", width = 900, height = 700)
dev.off()  # Close the PNG device

#############
# # Static: 'MS4A7','HLA-DQA1', 'DEFA4', 'BPI', 'CD1E', 'HLA-DRA', 'MS4A4E', 'SERPINB10', 'ELANE', 'CD1C'
# Top10: only HLA-DRAs express something but generally quite uniquely expressed in a certain cluster
FeaturePlot(object = integrated, features = c('MS4A7','HLA-DQA1', 'DEFA4', 'BPI', 'CD1E', 'HLA-DRA', 'MS4A4E', 'SERPINB10', 'ELANE', 'CD1C'))
dev.copy(png, file = "240701-top10-expression-static-1059x705.png", width = 1059, height = 705)
dev.off()  # Close the PNG device

DotPlot(object = integrated, features = c('MS4A7','HLA-DQA1', 'DEFA4', 'BPI', 'CD1E', 'HLA-DRA', 'MS4A4E', 'SERPINB10', 'ELANE', 'CD1C'), group.by = "Condition") +
  ggtitle("TOP10 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

dev.copy(png, file = "240701-top10-expression-static-dotplot-900x700.png", width = 900, height = 700)
dev.off()  # Close the PNG device

# Top20: C1QB, S100P
FeaturePlot(object = integrated, features = c('P2RY13','IL1RN','CD14','RETN','FGR','C1QB','MCEMP1','S100P','HLA-DRB5','AC016168.4'))
dev.copy(png, file = "240701-top20-expression-static-1059x705.png", width = 1059, height = 705)
dev.off()

DotPlot(object = integrated, features = c('P2RY13','IL1RN','CD14','RETN','FGR','C1QB','MCEMP1','S100P','HLA-DRB5','AC016168.4'), group.by = "Condition") +
  ggtitle("TOP20 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240701-top20-expression-static-dotplot-900x700.png", width = 900, height = 700)
dev.off()  # Close the PNG device

#############
# # Laminar: 'FAM167B', 'HSPA12B', 'CCL8', 'LYPD5', 'CAPN8', 'SERPINB4', 'KLF2', 'AC091563.1', 'KLHDC8A', 'LAMB3'
# Top10: FAM167B, HSPA12B, KLF2 expressed something but generally quite uniquely expressed in a 
# certain cluster that's different from Static
FeaturePlot(object = integrated, features = c('FAM167B', 'HSPA12B', 'CCL8', 'LYPD5', 'CAPN8', 'SERPINB4', 'KLF2', 'AC091563.1', 'KLHDC8A', 'LAMB3'))
dev.copy(png, file = "240701-top10-expression-laminar-1059x705.png", width = 1059, height = 705)
dev.off()  # Close the PNG device

DotPlot(object = integrated, features = c('FAM167B', 'HSPA12B', 'CCL8', 'LYPD5', 'CAPN8', 'SERPINB4', 'KLF2', 'AC091563.1', 'KLHDC8A', 'LAMB3'), group.by = "Condition") +
  ggtitle("TOP10 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

dev.copy(png, file = "240701-top10-expression-laminar-dotplot-900x700.png", width = 900, height = 700)
dev.off()  # Close the PNG device

# Top20: quite specific expression A2M, IGFBP5, COL15A1
FeaturePlot(object = integrated, features = c('FABP3','A2M','AL161629.1','AC110611.1','AL365214.2','IGFBP5','SLFN11','KLF4','COL15A1','AL136084.3'))
dev.copy(png, file = "240701-top20-expression-laminar-1059x705.png", width = 1059, height = 705)
dev.off() 

DotPlot(object = integrated, features = c('FABP3','A2M','AL161629.1','AC110611.1','AL365214.2','IGFBP5','SLFN11','KLF4','COL15A1','AL136084.3'), group.by = "Condition") +
  ggtitle("TOP20 GENES ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

dev.copy(png, file = "240701-top20-expression-laminar-dotplot-900x700.png", width = 900, height = 700)
dev.off()  # Close the PNG device


# # # # # # 
# JC CLUB #
# # # # # # 
# megakaryocytes?
# NO CD41, CD45, TPOR
# LOW TPO
FeaturePlot(object = integrated, features = c('CD36' ,'CXCR4', 'TPO', 'TRPC6')) 
dev.copy(png, file = "240701-TEST-1059x705.png", width = 1059, height = 705)
dev.off() 

# ERYTHROCYTES:
# https://www.researchgate.net/figure/Key-transcription-factors-involved-in-erythropoiesis-Transcription-factors-that-play_fig4_371366304
# NO CD71, CD235A, BND3
FeaturePlot(object = integrated, features = c('CD36' , 'BCL11A', 'SOX6','GYPA','SLC4A1',  'GATA1', 'KLF1', 'LDB1', 'MYB', 'TAL1')) 
dev.copy(png, file = "240701-CELL-TYPING2-1059x705.png", width = 1059, height = 705)
dev.off() 

DotPlot(object = integrated, features = c('CD36' , 'BCL11A', 'SOX6','GYPA','SLC4A1',  'GATA1', 'KLF1', 'LDB1', 'MYB', 'TAL1'), group.by = "Condition") +
  ggtitle("TRANSMEMBRANE TRANSPORTER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

dev.copy(png, file = "240701-CELL-TYPING2-dotplot-900x700.png", width = 900, height = 700)
dev.off()

# 

# TRANSMEMBRANE TRANSPORTERS
FeaturePlot(object = integrated, features = c('TRPC6', 'PIEZO1', 'PIEZO2', 'SLC4A11', 'KCNMA1', 'TRPV4', 'AQP1'))
dev.copy(png, file = "240701-TRANSMEMBRANE-TRANSPORTER-EXPRESSION-1059x705.png", width = 1059, height = 705)
dev.off() 

# LOW AQP5,8,9, LOW TRPM1
DotPlot(object = integrated, features = c('TRPC6', 'PIEZO1', 'PIEZO2', 'SLC4A11', 'KCNMA1', 'TRPV4', 'AQP1'), group.by = "Condition") +
  ggtitle("TRANSMEMBRANE TRANSPORTER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

dev.copy(png, file = "240701-TRANSMEMBRANE-TRANSPORTER-EXPRESSION-dotplot-900x700.png", width = 900, height = 700)
dev.off()  # Close the PNG device
# AQP <- related to vacuole formation during EHT, AQP5/8 = 0, otherwise all have low expression 
# ======================================================240701 ====================================================#


# # # # # # # # # # # # ##
# FROM LITERAUTRE REVIEW #
# # # # # # # # # # # # ##
# CD326 <- Erythroid, NOT DETECTED
# CD9 <- Megakaryocyte, universal expression
# CD19 <- leukocyte, few dots at the HP lineage
# no PG, HOXA, CD45, CD19, CD326, low CD7
FeaturePlot(object = integrated, features = c('YAP1', 'RUNX1', 'CD9'))





# https://cellregeneration.springeropen.com/articles/10.1186/s13619-024-00192-z
# # # GLUCOSE METABOLIC # # #
# HK1 has a strong binding affinity to mitochondria and functions as a housekeeping protein of glucose metabolism, 
# whereas HK2 binds less avidly to VDAC1 and alternates between cytoplasmic and mitochondrial-bound states in 
# response to environmental and metabolic stress
FeaturePlot(object = integrated, features = c('SLC2A1', 'SLC2A3','HK1','HK2'))
VlnPlot(integrated, features = c('SLC2A1', 'SLC2A3','HK1', 'HK2'))
DotPlot(object = integrated, features = c('SLC2A1', 'SLC2A3','HK1','HK2'), group.by = "Condition") +
  ggtitle("GLUCOSE METABOLIC") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(3,7)) +
  RotatedAxis()


# OXIDATIVE PHOSPHORYLATION - highly expressed across cluster
FeaturePlot(object = integrated, features = c('ATP5MC1', 'COX5A', 'CYC1'))
VlnPlot(integrated, features = c('ATP5MC1', 'COX5A', 'CYC1'))

DotPlot(object = integrated, features = c('ATP5MC1', 'COX5A', 'CYC1'), group.by = "Condition") +
  ggtitle("OXIDATIVE PHOSPHORYLATION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")



# GLYCOLYSIS - highly expressed across cluster
FeaturePlot(object = integrated, features = c('LDHA', 'TPI1', 'PKM'))
VlnPlot(integrated, features = c('LDHA', 'TPI1', 'PKM'))

DotPlot(object = integrated, features = c('LDHA', 'TPI1', 'PKM'), group.by = "Condition") +
  ggtitle("GLYCOLYSIS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()

DotPlot(object = integrated, features = c('PDHA1', 'PKM', 'LDHA', 'HIF1A','TPI1'), group.by = "Condition") +
  ggtitle("GLYCOLYSIS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# Glucose-oxidative-mesh
DotPlot(object = integrated, features = c('SLC2A1', 'SLC2A3','HK1','HK2','ATP5MC1', 'COX5A', 'CYC1', 'LDHA', 'TPI1', 'PKM'), group.by = "Condition") +
  ggtitle("METABOLIC GENES ACROSS CONDITIONS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# updated glycosis-oxidative shift
DotPlot(object = integrated, features = c('NRF2','MGST1','SOD2','ATP5MC1', 'COX5A', 'CYC1', 'PDHA1', 'PKM', 'LDHA', 'HIF1A','TPI1'), group.by = "Condition") +
  ggtitle("METABOLIC GENES ACROSS CONDITIONS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# DETOXIFICATION OF ROS - # nothing too interesting
FeaturePlot(object = integrated, features = c('BNIP3', 'PRDX2', 'PRDX3', 'PRDX4', 'P4HB', 'GPX7', 'ATOX1', 'NUDT2', 'CYCS', 'TXN2'))

DotPlot(object = integrated, features = c('BNIP3', 'PRDX2', 'PRDX3', 'PRDX4', 'P4HB', 'GPX7','ATOX1', 'NUDT2', 'CYCS', 'TXN2' ), group.by = "Condition") +
  ggtitle("DETOXIFICATION OF ROS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")


# FATTY ACID METABOLIC
# universally expressed across cluster, didn't save the img
FeaturePlot(object = integrated, features = c('SDHC', 'CA2', 'ADSL', 'ECHS1', 'DLD', 'ECI2'))
# S100A10 exp (high to low): ENdo > lineage > bottleneck, saved img
FeaturePlot(object = integrated, features = c('ODC1', 'IDH1', 'S100A10', 'CCDC58'))

DotPlot(object = integrated, features = c('SDHC', 'CA2', 'ADSL', 'ECHS1', 'DLD', 'ECI2','ODC1', 'IDH1', 'S100A10', 'CCDC58' ), group.by = "Condition") +
  ggtitle("FATTY ACID METABOLIC") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# GLUTAMINE METABOLIC
# 'ASL', 'GOT1' GOT REMOVED DUE TO NON-SPECIFIC EXPRESSION
# FAH gpt upregulated in static and pulsatile
# PYCR1 specically up-regulated in all condition but prestim
FeaturePlot(object = integrated, features = c('GAD1', 'FAH', 'PYCR1', 'NIT2'))
# doesn't show much difference across cluster, img not saved
FeaturePlot(object = integrated, features = c('GOT2', 'ADSL', 'GPT2', 'GLUD1'))

DotPlot(object = integrated, features = c('GAD1', 'FAH', 'PYCR1', 'NIT2', 'GOT2', 'ADSL', 'GPT2', 'GLUD1' ), group.by = "Condition") +
  ggtitle("GLUTAMINE METABOLIC") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# ================================================================================ #
# # # # # # # #
# CELL TYPING #
# # # # # # # #
# ENDOTHELIAL CELL MARKER - looks fine
# https://www.proteinatlas.org/ENSG00000187513-GJA4/tissue+cell+type
FeaturePlot(object = integrated, features = c('DLL4','GJA5', 'EFNB2', 'HEY1', 'NOTCH4', 'SOX17', 'EPHB4', 'NR2F2', 'FLT4'))

DotPlot(object = integrated, features = c('DLL4','GJA5', 'EFNB2', 'HEY1', 'NOTCH4', 'SOX17', 'EPHB4', 'NR2F2', 'FLT4'), group.by = "Condition") +
  ggtitle("ENDOTHELIAL CELL MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")


# ARTERIAL: 'DLL4', 'GJA4', 'GJA5', 'EFNB2', 'SOX17', 'HEY1'
# HEY2 is more smooth muscle related
# GJA4 leaked passed bottleneck
# FeaturePlot(object = integrated, features = c('DLL4', 'GJA4', 'GJA5', 'EFNB2', 'SOX17', 'HEY1'))
FeaturePlot(object = integrated, features = c('DLL4','GJA4', 'GJA5', 'EFNB2', 'HEY1', 'NOTCH4', 'SOX17'))
DotPlot(object = integrated, features = c('DLL4','GJA5', 'EFNB2', 'HEY1', 'NOTCH4', 'SOX17'), group.by = "Condition") +
  ggtitle("ENDOTHELIAL CELL (ARTERIAL) MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# Venous: 'NR2F2', 'EPHB4', 'SMARCA4', 'FLT4'
# a bit leakage of EMCN towards the bottleneck and PASSED the bottleneck so removed
# FeaturePlot(object = integrated, features = c('NR2F2', 'EPHB4', 'FLT4'))
FeaturePlot(object = integrated, features = c('EPHB4', 'NR2F2', 'FLT4'))
DotPlot(object = integrated, features = c('EPHB4', 'NR2F2', 'FLT4'), group.by = "Condition") +
  ggtitle("ENDOTHELIAL CELL (VENOUS) MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# EHT marker: 'RUNX1', 'SOX17', 'GATA2', 'MEIS1', 'GFI1', 'CD44' (CD41, CD43, CD45,  'CD31', 'VEC', 'FLK-1' not founc)
# GFI1/ GFI1B: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9136496/
# MEIS1 (earliest EHT marker + prior to RUNX1): https://pubmed.ncbi.nlm.nih.gov/37500618/
# EHT markers: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6490701/
# I don't understand what's causing no/low SOX17 expression pass the endothelial cells
# not found: CD41, CD31, VEC, FLK-1 
FeaturePlot(object = integrated, features = c('RUNX1', 'SOX17', 'GATA2', 'MEIS1', 'GFI1', 'GFI1B'))
DotPlot(object = integrated, features = c('RUNX1', 'SOX17', 'GATA2', 'MEIS1', 'GFI1', 'GFI1B'), group.by = "Condition") +
  ggtitle("EHT MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")


# Haematopoietic lineage markers: 'KLF1', 'CD14', 'ITGA2B', 'MPO', 'CSF2RA'
# 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'CD1A' 
# klf1, ITGA2B: erythoid
FeaturePlot(object = integrated, features = c('KLF1','ITGA2B', 'CD14',  'MPO', 'CSF2RA'))
DotPlot(object = integrated, features = c('KLF1','ITGA2B', 'CD14',  'MPO', 'CSF2RA'), group.by = "Condition") +
  ggtitle("HAEMATOPOIETIC LINEAGE MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")



# lymphoid: https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.b.21614
# T, B, NK CELLS <- CD3, CD19, CD56, HLA-DR, CD34 
# CD3, CD56 NOT FOUND
# cd19: all cells have same value
FeaturePlot(object = integrated, features = c('HLA-DRA', 'HLA-DRB1', 'HLA-DRB5'))
# all at static
DotPlot(object = integrated, features = c('HLA-DRA', 'HLA-DRB1', 'HLA-DRB5'), group.by = "Condition") +
  ggtitle("HAEMATOPOIETIC LINEAGE MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# myeloid: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4742930/; https://www.biocompare.com/Editorial-Articles/612866-A-Guide-to-Myeloid-Cell-Markers/
# not found: CD45, CD206, CD169, CD11c, CD123, CD16, CD15, Ly6C, CD1a, CD141, CD11b, CD64, CD1c, CD103, CD11b
# NON-SPECIFIC PASS BOTTLENECK ARE EXCLUDED: 'ANGPT2','ARG1', 'CCL2','CSF1','ANPEP', 'ADGRE1','CD24', IL3RA
# LOW EXPRESSION: CD1A
# IL3RA/ CD123: responsible for the proliferation, survival, and differentiation of hematopoietic cells.
# MPO: peroxidase enzyme, != myloid
# CSF2RA: granulocytes and macrophages
FeaturePlot(object = integrated, features = c( 'CX3CR1', 'CSF2RA', 'CD14', 'CCR2', 'MPO'))
DotPlot(object = integrated, features = c('CX3CR1', 'CSF2RA', 'CD14', 'CCR2', 'MPO'), group.by = "Condition") +
  ggtitle("MYLEOID MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# HAEMATOPOEITIC LINEAGE MASH
FeaturePlot(object = integrated, features = c( 'CX3CR1', 'CSF2RA', 'CD14', 'CCR2', 'CD1A', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'KLF1', 'ITGA2B', 'MPO'))
DotPlot(object = integrated, features = c('CX3CR1', 'CSF2RA', 'CD14', 'CCR2', 'CD1A', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'KLF1', 'ITGA2B', 'MPO'), group.by = "Condition") +
  ggtitle("HAEMATOPOEITIC LINEAGE") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")


# =============================== TEST ===================================#
# AQP5 = 0
# AQP9 = 0/ super low
FeaturePlot(object = integrated, features = c('AQP1', 'AQP3'))

DotPlot(object = integrated, features = c('AQP1', 'AQP3'), group.by = "Condition") +
  ggtitle("AQUAPORIN EXPRESSION ACROSS CONDITIONS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

# Define a named vector of colours, same colours as thowed on the heatmap
condition.colour <- c('Prestimulation' = '#F8766D',
                      'Static' = '#C2D989',
                      'Laminar' = '#C476FF')
# co-localisation by conditions
FeatureScatter(integrated, 
               'AQP1', 'AQP3', 
               group.by = "Condition", 
               pt.size = 1.5, 
               cols = condition.colour)
# co-localisation in one pot
FeatureScatter(integrated, 
               'AQP1', 'AQP3',
               pt.size = 3,
               cols = condition.colour) +
  facet_wrap(Idents(integrated))

DotPlot(object = integrated, features = c('YAP1', 'LATS1', 'RUNX1'), group.by = "Condition") +
  ggtitle("EHT") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis() +
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
