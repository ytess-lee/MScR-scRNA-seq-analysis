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
integrated$Condition <- factor(integrated$Condition, levels = c("Prestimulation", "Static", "Laminar"))

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


# # # # # # 
# JC CLUB #
# # # # # # 
# megakaryocytes?
# NO CD41, CD45, TPOR
# LOW TPO
FeaturePlot(object = integrated, features = c('CD36' ,'CXCR4', 'TPO', 'TRPC6')) 
dev.copy(png, file = "240630-TEST-1059x705.png", width = 1059, height = 705)
dev.off() 

# ERYTHROCYTES:
# https://www.researchgate.net/figure/Key-transcription-factors-involved-in-erythropoiesis-Transcription-factors-that-play_fig4_371366304
# NO CD71, CD235A, BND3
FeaturePlot(object = integrated, features = c('CD36' , 'BCL11A', 'SOX6','GYPA','SLC4A1',  'GATA1', 'KLF1', 'LDB1', 'MYB', 'TAL1')) 
dev.copy(png, file = "240630-CELL-TYPING2-1059x705.png", width = 1059, height = 705)
dev.off() 

l1 <- DotPlot(object = integrated, features = c('CD36' , 'BCL11A', 'SOX6','GYPA','SLC4A1',  'GATA1', 'KLF1', 'LDB1', 'MYB', 'TAL1'), group.by = "Condition") +
  ggtitle("TRANSMEMBRANE TRANSPORTER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

dev.copy(png, file = "240630-CELL-TYPING2-dotplot-900x700.png", width = 900, height = 700)
dev.off()

# 

# TRANSMEMBRANE TRANSPORTERS
FeaturePlot(object = integrated, features = c('LRRC8A','SLC4A11', 'KCNMA1', 'TRPV4','TRPC6', 'PIEZO2','PIEZO1', 'AQP1'))
dev.copy(png, file = "240707-TRANSMEMBRANE-TRANSPORTER-EXPRESSION-1059x705.png", width = 1059, height = 705)
dev.off() 

# LOW AQP5,8,9, LOW TRPM1 (removed)
DotPlot(object = integrated, features = c( 'LRRC8A','SLC4A11', 'KCNMA1', 'TRPV4','TRPC6', 'PIEZO2','PIEZO1', 'AQP1'), group.by = "Condition") +
  ggtitle("TRANSMEMBRANE TRANSPORTER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

dev.copy(png, file = "240707-TRANSMEMBRANE-TRANSPORTER-EXPRESSION-dotplot-800x700.png", width = 800, height = 700)
dev.off() 
