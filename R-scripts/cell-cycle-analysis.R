# this script is for cell cycle analysis
# first group clusters by conditions (Prestim, static, laminar)
# check cell cycle gene expression

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
# https://cyclebase.org/Advanced%20Search?type=9606
# G1 peak time: 'CCNE1','NXF1','TRMT2A','DNAJB6','INTS8'
# G1-S1 PHASE: 'CCNE2','CDK2'
# S PEAK TIME: 'PCNA','RRM2', 'TOPBP1'+ 'CCND3', 'CDK4'
# S-G2 TRANSITION: 'CCNA2'
# G2: 'CDK1'
# G2-M TRANSITION: 'CCNB1'
# M PHASE: 'PLK1','NUF2','RAD21','CKAP5', 'CDC27'


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


integrated <- readRDS("240719-integrated-with-cc-data.RDS") 
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.3)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis

temp_idents <- integrated@active.ident
Idents(integrated) <- "Condition" 
integrated$Condition <- factor(integrated$Condition, levels = c("Prestimulation", "Static", "Laminar"))

d13 <- DimPlot(integrated, reduction = "umap", label.size = 5)

# G1 peak time: 'CCNE1','NXF1','TRMT2A','DNAJB6','INTS8'
# G1-S1 PHASE: 'CCNE2','CDK2' 
# low CCNE1, NXF1, TRMT2A
DotPlot(integrated, features = c('CCNE1','NXF1','TRMT2A','DNAJB6','INTS8',
                                 'CCNE2','CDK2' )) +
  ggtitle("G1 + G1-S1 PHASE MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Features") +
  ylab("Conditions")

Idents(integrated) <- temp_idents
Idents(integrated) <- 'Condition'
integrated$Condition <- factor(integrated$Condition, levels = c("Prestimulation", "Static", "Laminar"))
DotPlot(integrated, features = c('CDK4','CDK6', # G1
                                 'CCNE2','CDK2', # G1/S checkpoint
                                 'CCND3', # S
                                 'CCNA2','CDK1', # S/G2 transition
                                 'CCNB1' )) + # G2M transition
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab("") +
  ylab("")
dev.copy(png, file = 'output-dissertation-figures_240720/figure2/240804-integrated-cc-markers-condition-480x200.png', width = 480, height = 200)
dev.off()

# CELL CYCLE ANALYSIS (PREV_SCRIPT + SEURAT VIN)
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
# Seurat.forCC <- CellCycleScoring(integrated, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE) # if run without change to different scaling method, warning: the following features are not present in the object: MLF1IP/ FAM64A, HN1, not searching for symbol synonyms

Seurat.forCC <- NormalizeData(integrated)
Seurat.forCC <- FindVariableFeatures(Seurat.forCC, selection.method = "vst")
Seurat.forCC <- ScaleData(Seurat.forCC, features = rownames(Seurat.forCC))

# saveRDS(Seurat.forCC, file = '240719-seurat-obj-from-integrated-for-cell-cycle-analysis.RDS')
Seurat.forCC <- readRDS(file = '240719-seurat-obj-from-integrated-for-cell-cycle-analysis.RDS')
Idents(Seurat.forCC)
Seurat.forCC <- CellCycleScoring(Seurat.forCC, s.features = s_genes, 
                                 g2m.features = g2m_genes, set.ident = TRUE)
# wait running this line gives the same warning message as running it stright with SCT
# Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
# Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms

DimPlot(Seurat.forCC, split.by = "Condition", ncol = 2, pt.size = 0.1) +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) 

dev.copy(png, file = '240719-integrated-cc-analysis-condition-856x676.png', width = 856, height = 676)
dev.off()

DimPlot(subset(Seurat.forCC, idents = "seurat_clusters" ))

# Move cell cycle data back into Seurat object
integrated$CellCycle <- Idents(Seurat.forCC)
integrated$CellCycle <- factor(integrated$CellCycle, levels = c("G1", "S", "G2M"))
integrated$S.Score <- Seurat.forCC@meta.data$S.Score
integrated$G2M.Score <- Seurat.forCC@meta.data$G2M.Score

# create a table with the cell cycle score and the library ID
tb <- table(integrated$CellCycle, integrated$Condition)

# calculate the percentage of the various cell cycle score
tb <- apply(tb, 1, function(x){x/colSums(tb)*100})

# create a data frame to use for the plotting in GG plot
perc.cc <- data.frame(tb) %>% 
  rownames_to_column(var = "Cluster") %>% 
  pivot_longer(-Cluster, names_to = "Stage", values_to = "Percentage")

perc.cc$Stage <- factor(perc.cc$Stage, levels = c("G1", "S", "G2M"))

DimPlot(integrated, split.by = "CellCycle")+
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) 
dev.copy(png, file = '240719-integrated-cc-analysis-10clusters-944x463.png', width = 944, height = 463)
dev.off()


DimPlot(Seurat.forCC, split.by = "Condition") +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) 
dev.copy(png, file = '240719-integrated-cc-analysis-condition-944x463.png', width = 944, height = 463)
dev.off()

saveRDS(integrated, file = '240719-integrated-with-cc-data.RDS')

ggplot(perc.cc, aes(x = Stage, y = Percentage,)) +
  geom_col(aes(fill = Cluster), position="dodge") +
  theme(text = element_text(size=13),
        axis.title = element_text(size = 16))
dev.copy(png, file = '240719-integrated-cc-analysis-condition-barplot-944x463.png', width = 944, height = 463)
dev.off()

temp_cc_cluster <- perc.cc$Cluster
perc.cc$Cluster <- factor(perc.cc$Cluster, levels = c("Prestimulation", "Static", "Laminar"))

ggplot(perc.cc, aes(x = Cluster, y = Percentage)) +
  geom_col(aes(fill = Stage)) +
  theme(text = element_text(size = 16),
    axis.title = element_text(size = 16)) + RotatedAxis() + ggtitle(NULL) + xlab('') + ylab('')
dev.copy(png, file = 'output-dissertation-figures_240720/figure2/240804-integrated-cc-analysis-condition-barplot(in-one)-300x540.png', width = 300, height = 540)
dev.off()

# theme(
#   text = element_text(family = "Arial", size = 12),  # Change 'Arial' to desired font type, and 12 to desired font size
#   axis.title = element_text(size = 16),  # Change 16 to desired font size for axis titles
#   axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text if necessary
  

# colocalisation by clusters
FeatureScatter(integrated, "RUNX1", "S.Score") +
  facet_wrap(Idents(integrated))
dev.copy(png, file = '240719-integrated-cc-analysis-RUNX1-S.Score-colocalisation-944x463.png.png', width = 944, height = 463)
dev.off()

FeatureScatter(integrated, "RUNX1", "G2M.Score") +
  facet_wrap(Idents(integrated))
dev.copy(png, file = '240719-integrated-cc-analysis-RUNX1-G2M.Score-colocalisation-944x463.png.png', width = 944, height = 463)
dev.off()

FeatureScatter(integrated, "RUNX1", "MYBL2")+
  facet_wrap(Idents(integrated)) 
dev.copy(png, file = '240719-integrated-cc-analysis-RUNX1-MYBL2-colocalisation-944x463.png.png', width = 944, height = 463)
dev.off()

#remove the object that we used to calculate the cell cycle score
rm(Seurat.forCC)

# Colocalisation with by condition
temp_idents <- integrated@active.ident
Idents(integrated) <- "Condition" 

FeatureScatter(integrated, "RUNX1", "S.Score") +
  facet_wrap(Idents(integrated))
dev.copy(png, file = '240719-integrated-cc-analysis-RUNX1-S.Score-colocalisation(by-condition)-944x463.png.png', width = 944, height = 463)
dev.off()

FeatureScatter(integrated, "RUNX1", "G2M.Score") +
  facet_wrap(Idents(integrated))
dev.copy(png, file = '240719-integrated-cc-analysis-RUNX1-G2M.Score-colocalisation(by-condition)-944x463.png.png', width = 944, height = 463)
dev.off()

FeatureScatter(integrated, "RUNX1", "MYBL2")+
  facet_wrap(Idents(integrated)) 
dev.copy(png, file = '240719-integrated-cc-analysis-RUNX1-MYBL2-colocalisation(by-condition)-944x463.png.png', width = 944, height = 463)
dev.off()


rm(Seurat.forcc)
Idents(integrated) <- "CellCycle" 
RidgePlot(integrated, features = c('CDK4','CDK6','CCNE2',
                                   'CDK2','CCND3', 'CCNA2',
                                   'CDK1', 'CCNB1'), ncol = 2, log = T, group.by = 'Condition')
dev.copy(png, file = '240719-integrated-cc-analysis-CC-genes(by-Condition)+log-1500x2000.png.png', width = 1500, height = 2000)
dev.off()

































# let's stick to cyberbase3.0 first, then proceeded to this paper for further analysis
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-017-00039-z/MediaObjects/41467_2017_39_MOESM1_ESM.pdf
# G1: 'DTL','PTTG1','CDKN3','ZNF367','SLBP','MCM6','HSPA8','CDCA7','SKP2','ANTXR1','IVNS1ABP','DYNLL1','GRPEL1','ZRANB2','OPN3', 'MSL1', 'AOC2'

FeaturePlot(integrated, features = c('DTL','PTTG1','CDKN3','ZNF367'), split.by = 'Condition')
VlnPlot(integrated, features = c('DTL','PTTG1','CDKN3','ZNF367'))

DotPlot(integrated, features = c('DTL','PTTG1','CDKN3','ZNF367','SLBP',
                                 'MCM6','HSPA8','CDCA7','SKP2','ANTXR1',
                                 'IVNS1ABP','DYNLL1','GRPEL1','ZRANB2','OPN3', 'MSL1', 'AOC2')) +
  ggtitle("G1 MARKERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Features") +
  ylab("Conditions")

# S:

# G2:

# M: 