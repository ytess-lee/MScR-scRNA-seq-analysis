# use res=0.3 for clustering
# then subset the top right cluster out for further analysis 
# |

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


integrated <- readRDS("RDS/240615_pleasework_Mechanobiology_Harmonised_2reps.RDS") 
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.3)
integrated <- RunUMAP(object = integrated, dims = 1:13) # it flips on the x-axis


# confirm that the subset works
integrated_umap <- DimPlot(integrated, reduction = "umap", label = T)
cluster_main_umap <- DimPlot(cluster_main, reduction = 'umap', label = T)
integrated_umap|cluster_main_umap
dev.copy(png, file = "240713-all-clusters(9)-subset-cluster_main(0-7)-1095x700.png", width = 1095, height = 700)
dev.off()

# merged UAMP
# umap <- DimPlot(integrated, reduction = "umap", label = T)
# umap
# 
# cluster_condition <- DimPlot(integrated, reduction = 'umap', split.by = 'Condition')
# cluster_condition
# 
# # temporarly store the cluster idents to be able to re-use them later
temp_idents <- integrated@active.ident
# 
# # assign the culture condition to the active identitites for the DEG analysis
Idents(integrated) <- "Condition" 
# 
# # Reorder the levels of the condition factor
integrated$Condition <- factor(integrated$Condition, levels = c("Prestimulation", "Static", "Laminar"))
# 
# by_condition <- DimPlot(integrated, reduction = 'umap', split.by = 'Condition')
# dev.copy(png, file = "240709-Condition(res0.3)-937x708.png", width = 937, height = 708)
# dev.off()
# 
# cluster_condition|by_condition
# dev.copy(png, file = "240707-split-by-day-con(res0.3)-2000x1000.png", width = 2000, height = 1000)
# dev.off()

# once resolution is determined, subset the top right clusters out (#8, 9)

# # ======================== # SUBSET CLUSTERS # ======================== # 
Idents(integrated) <- "SCT_snn_res.0.3"
cluster_endo <- subset(integrated, idents= c('0','1', '2','3','4','5'))
cluster8_9 <- subset(integrated, idents= c("8", "9"))
# # cluster6789 <- subset(integrated, idents= c("6", "7", "8", "9")) # no diff between just 8+9/ 6+8/ 6-8
# # cluster6_8_blood <- subset(integrated, idents = c('6', '8'))
# cluster_main <- subset(integrated, idents = c('0','1', '2','3','4','5','6','7'))
# cluster_7 <- subset(integrated, idents = '7')
# cluster6_7 <- subset(integrated, idents = c('6','7'))

# save RDS
# saveRDS(cluster_main, file = '240713-cluster0-7-subset.RDS')
# saveRDS(cluster8_9, file = '240709-cluster-8+9-subset.RDS')
# saveRDS(cluster6_8_blood, file = '240709-cluster-6+8-subset.RDS')
# saveRDS(cluster_7, file = '240714-cluster7-subset.RDS')
# saveRDS(cluster6_7, file = '240714-cluster6+7(EHT+Blood)-subset.RDS')
saveRDS(cluster_endo, file = '240718-cluster0-5(Endothelium)-subset.RDS')
# ======================== # SUBSET CLUSTER endo # ======================== # 
# once found the top10, GO analysis to do the cell typing
cluster_endo <- readRDS(file = '240718-cluster0-5(Endothelium)-subset.RDS')

de_endo <- FindAllMarkers(cluster_endo, only.pos = T) # for gene set ENRICHMENT analysis
filtered_markers <- de_endo %>% filter(p_val_adj < 0.05)

# write.csv(filtered_markers, "240718_filtered(p)_markers_new-data-clusters_endo(0-5)-subset.csv", row.names = FALSE)
filtered_markers <- read.csv(file = '240718_filtered(p)_markers_new-data-clusters_endo(0-5)-subset-AUTOMATED.csv')

top10_main <- filtered_markers  %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 


DoHeatmap(cluster_endo, features = top10_main$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240718-CLUSTER_main-subset-DEGs-1059x705.png", width = 1059, height = 705)
dev.off()
write.csv(top10_main, "240718_filtered(p)_markers_new-data-cluster0-5-Top10(omitted-removed).csv", row.names = FALSE)


DimPlot(cluster_main, reduction = 'umap')

DimPlot(integrated, reduction = "umap")
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = TOP10 = = = = = = = == = = = = = 
# cluster0+1: "CARTPT", "ACSM3", "AC007422.1", "MYB", "ACSM1", "LNCAROD", "SPINK2", "MIR924HG", "MKI67", "PIF1", "IGFBP5", "FAM167B", "NKAIN2", "CYSLTR1", "NCALD", "A2M", "KLF2", "ABCC3", "SGCZ", "THSD7A",
# cluster2+3: "FCGR2A", "GFI1B", "UBE2C", "NKG7", "CDC20", "RAD51AP1", "CDC45", "HIST1H1D", "GINS2", "BIRC5", "SST", "KRT7", "REN", "NPPB", "ACTG2", "ACTC1", "MYL7", "PRR9", "COL1A1", "APOA1", 
# cluster4+5: "AC114757.1", "CHRM3", "AC084740.1", "PRND", "TMEM132B", "AL356124.1", "SNTB1", "AC093817.2", "LRRC36", "KCNIP1", "FABP3", "MMP1", "KLF2", "FAM167B", "IL32", "NQO1", "SLC9A3R2", "C12orf57", "S100A6", "CCL2"
# cluster6+7: "PINCR", "AL023574.1", "AC044893.1", "HS3ST4", "PURPL", "PLD5", "MIR3659HG", "LINC01707", "NELL1", "MS4A6A", "ADAMTS18", "DSCAM", "ADAMTSL1", "CXCL8", "PTX3", "TNFRSF11B", "ANKRD1", "EDN1", "CXCL12", "CXCL1"

FeaturePlot(cluster_main, features = c("CARTPT", "ACSM3", "AC007422.1", "MYB", "ACSM1", "LNCAROD", "SPINK2", "MIR924HG", "MKI67", "PIF1", "IGFBP5", "FAM167B", "NKAIN2", "CYSLTR1", "NCALD", "A2M", "KLF2", "ABCC3", "SGCZ", "THSD7A"))
dev.copy(png, file = '240713-Top10-cluster0+1-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(cluster_main, features = c("FCGR2A", "GFI1B", "UBE2C", "NKG7", "CDC20", "RAD51AP1", "CDC45", "HIST1H1D", "GINS2", "BIRC5", "SST", "KRT7", "REN", "NPPB", "ACTG2", "ACTC1", "MYL7", "PRR9", "COL1A1", "APOA1"))
dev.copy(png, file = '240713-Top10-cluster2+3-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(cluster_main, features = c("AC114757.1", "CHRM3", "AC084740.1", "PRND", "TMEM132B", "AL356124.1", "SNTB1", "AC093817.2", "LRRC36", "KCNIP1", "FABP3", "MMP1", "KLF2", "FAM167B", "IL32", "NQO1", "SLC9A3R2", "C12orf57", "S100A6", "CCL2"))
dev.copy(png, file = '240713-Top10-cluster4+5-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(cluster_main, features = c("PINCR", "AL023574.1", "AC044893.1", "HS3ST4", "PURPL", "PLD5", "MIR3659HG", "LINC01707", "NELL1", "MS4A6A", "ADAMTS18", "DSCAM", "ADAMTSL1", "CXCL8", "PTX3", "TNFRSF11B", "ANKRD1", "EDN1", "CXCL12", "CXCL1"))
dev.copy(png, file = '240713-Top10-cluster6+7-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

# Group 0: Proliferating metabolic cells (KEGG - CELL CYCLE)
# Group 1: Endothelial cells (KEGG - Oxidative phosphorylation)
# Group 2: Proliferating cells. (KEGG - Oxidative phosphorylation)
# Group 3: Epithelial and muscle cells (kidney/heart). (KEGG - Oxidative phosphorylation)
# Group 4: Neural cells (KEGG - immune cell differentiation)
# Group 5: Stromal cells (KEGG - Oxidative phosphorylation)
# Group 6: Haematopoietic lineage
# Group 7: EHT
# Group 8: Megakaryocyte
# Group 9: Ribosome biogenesis


# AF blood naive 1:'SOX4','LM04','ID2','MEIS2','TSC22D1','GATA2','SPI1','JUNB','BHLHE40','MITF','CEBPB','HOXB2','KDMSB','TGIF1','AHR','MSRB2','GLMP'
# Naive 2 (-naive 1): 'ID4','HES1','FOXP1','HMGB2','MITF','GFI1')
# I don't think it gives out much info about the endothelial population
FeaturePlot(integrated, features = c('SOX4','ID2','MEIS2','TSC22D1','GATA2','SPI1','JUNB','BHLHE40','MITF','CEBPB','HOXB2','TGIF1','AHR','MSRB2','GLMP', 'ID4','HES1','FOXP1','HMGB2','MITF','GFI1'))
dev.copy(png, file = '240713-Naive-markers-fromAF2020-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

# ery: 'KLF1','MYC','HMGB1','HMGB2','TFAM','ZFP35L2','PA2G4','NOLC1','YBX1','PTTG1'
# gran: 'CEBPD','SPI1','
FeaturePlot(integrated, features = c('KLF1','MYC','HMGB1','HMGB2','TFAM','ZFP35L2','PA2G4','NOLC1','YBX1','PTTG1'))
dev.copy(png, file = '240713-ERY-markers-fromAF2020-2000x2000.png', width = 2000, heigh = 2000)
dev.off()


# ======================== # SUBSET CLUSTER 8+9 # ======================== # 
cluster8_9 <- readRDS(file = 'RDS/240709-cluster-8+9-subset.RDS')

de_markers_condition <- FindAllMarkers(cluster8_9, only.pos = T) # for gene set ENRICHMENT analysis
filtered_markers <- de_markers_condition %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, "240709_filtered(p)_markers_new-data-clusters8+9-subset.csv", row.names = FALSE)
#my.markers <- read.csv(file = '240709_filtered(p)_markers_new-data-clusters8+9-subset.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) # top 10/20 just for the heatmap 

# write.csv(top10, file ="240519-top10-avg_log2FC.csv")
sum(top10$gene %in% rownames(GetAssayData(cluster8_9, layer = 'scale.data')))
DoHeatmap(cluster8_9, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240709-CLUSTER8+9-subset-DEGs-1059x705.png", width = 1059, height = 705)
dev.off()



# to see if i could further differentiate the clusters from 8+9
# apparently you can
# tried with high res but the heatmap doesnt show representatitive markers between 0 and 1 (0,1,2) so stick to 0+1 (res=0.1)
DimPlot(cluster8_9)

cluster8_9 <- FindClusters(cluster8_9, dims= 1:15, resolution = 0.1)
DimPlot(cluster8_9, reduction = "umap") # separate into 2 clusters
cluster8_9 <- RunUMAP(object = cluster8_9, dims = 1:15) # it flips on the x-axis
DimPlot(cluster8_9, reduction = 'umap', label=T)

# mega: 'GP9', 'PF4', 'GATA1', 'TAL1', 'FLI1'
# ery: 'GYPA', 'KLF1','MYC', 
# gran:''AZU1', CEBPD', 'CEBPB', 'CEBPA', 'CEBPE'
feature <- FeaturePlot(cluster8_9, features = c( 'GP9', 'PF4', 'GATA1', 'TAL1', 'FLI1', 'GYPA', 'KLF1','MYC', 'HBZ', 'HBM'))
feature
dev.copy(png, file = "240714-CLUSTER8+9-subset-ery+mega-1059x800.png", width = 1059, height = 800)
dev.off()

dim <- DotPlot(cluster8_9, features = c('GP9', 'PF4', 'GATA1', 'TAL1', 'FLI1', 'GYPA', 'KLF1','MYC', 'HBZ', 'HBM'), group.by = "seurat_clusters") +
  ggtitle("MEGAKARYOCTYE-/ ERYTHROID-COMMITTED MARKER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Features") +
  ylab("Clusters")
dim
dev.copy(png, file = "240714-CLUSTER8+9-subset-ery+mega-dotplot-900x700.png", width = 500, height = 700)
dev.off()

de89 <- FindAllMarkers(cluster8_9, only.pos = T) # for gene set ENRICHMENT analysis
filtered_markers <- de89 %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, "240714_filtered(p)_markers_new-data-clusters8+9(res0.1+2-sub-clusters)-subset.csv", row.names = FALSE)
filtered_markers <- read.csv(file= 'output-subset-1-7, 8+9/8+9(static-outliner)/csv/240714_filtered(p)_markers_new-data-clusters8+9(res0.1+2-sub-clusters)-subset-AUTOMATED.csv')
top10_89 <- filtered_markers  %>% 
  group_by(cluster) %>% 
  top_n(20, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 

DoHeatmap(cluster8_9, features = top10_89$gene) + theme(axis.text.y = element_text(size = 10))
dev.copy(png, file = "240709-CLUSTER8+9(0+1)-subset-DEGs-1059x705.png", width = 1059, height = 705)
dev.off()

# shinygo
# cluster 0: Ery + Mega
# cluster 1: dunno, highly expressing gene related to Hb assembly
# there's a population within cluster 1 that doesn't express top10 markers from both clusters
# "AL139815.1", "AHSP", "AC100801.1", "AC026167.1", "AC024267.1", "AC011246.1", "AC008415.1", "AC004053.1"
      # no actually couldn't tell by eye

excluded <- FeaturePlot(cluster8_9, features = c("AL139815.1", "AHSP", "AC100801.1", "AC026167.1", "AC024267.1", "AC011246.1", "AC008415.1", "AC004053.1"))
excluded
dev.copy(png, file = "240714-CLUSTER8+9(0+1)-subset-subgroup-within-cluster1-1059x705.png", width = 1059, height = 705)
dev.off()



# ======================== # SUBSET CLUSTER 7 # ======================== # 
cluster_7 <- readRDS(file = '240714-cluster7-subset.RDS')
cluster_7 <- FindClusters(cluster_7, dims= 1:15, resolution = 0.1)
c7 <- DimPlot(cluster_7, reduction = "umap") 
original <- DimPlot(integrated, reduction = 'umap')
c7|original

# de_markers_7 <- FindAllMarkers(cluster_7, only.pos = T) # for gene set ENRICHMENT analysis
# filtered_markers <- de_markers_7 %>% filter(p_val_adj < 0.05)
# write.csv(filtered_markers, "240714_filtered(p)_markers_new-data-clusters7-subset.csv", row.names = FALSE)

filtered_markers <- read.csv(file = '240714_filtered(p)_markers_new-data-clusters7-subset-AUTOMATED.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) # top 10/20 just for the heatmap 

# write.csv(top10, file ="240519-top10-avg_log2FC.csv")

DoHeatmap(cluster_7, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240714-CLUSTER7-subset-DEGs-1059x705.png", width = 1059, height = 705)
dev.off()

#bottleneck
FeaturePlot(cluster_7,features = c('GATA2', 'ITGA2B', 'CD44' ,'RUNX1' ,'CD82','DOK2', 'LGALS9', 'PTK2B','MCTP1', 'CORO1A', 'SPI1' ))
#blood
FeaturePlot(cluster_7, features = c('NFE2', 'PTK2B', 'CORO1A', 'SPI1' , 'MYB', 'KLF1', 'GATA1', 'GFI1','TRPV4' ))  
#endo
FeaturePlot(cluster_7, features = c('TEK','PCDH12', 'RAMP2', 'HEY1','HOPX', 'MECOM','DLL4', 'KRT18', 'NR2F2', 'SOX17'))

c7
# cluster 0 (MID): 'MECOM','HOPX','KRT18','RAMP2'
# CLUSTER 1 (RIGHT): 'CD44','RUNX1','CD82','SPI1','NFE2','PTK2B','CORO1A','MYB'
# CLUSTER 2 (LEFT): 'TEK','PCDH12','SOX17','HEY1','DLL4'
f <- FeaturePlot(cluster_7, features = c('TEK','PCDH12', 'RAMP2', 'HEY1','HOPX', 'MECOM','DLL4', 'KRT18', 'NR2F2', 'SOX17',
                                    'GATA2', 'ITGA2B', 'CD44' ,'RUNX1' ,'CD82','DOK2', 'LGALS9', 'MCTP1', 'SPI1',
                                    'NFE2', 'PTK2B', 'CORO1A' , 'MYB', 'KLF1', 'GATA1', 'GFI1','TRPV4' ))
dev.copy(png, file = "240714-CLUSTER7-subset-3000x5000.png", width = 3000, height = 5000)
dev.off()
c7

DotPlot(cluster_7, features = c('TEK','PCDH12', 'RAMP2', 'HEY1','HOPX', 'MECOM','DLL4', 'KRT18', 'NR2F2', 'SOX17',
                                    'GATA2', 'ITGA2B', 'CD44' ,'RUNX1' ,'CD82','DOK2', 'LGALS9','MCTP1', 'SPI1',
                                    'NFE2', 'PTK2B', 'CORO1A', 'MYB', 'KLF1', 'GATA1', 'GFI1','TRPV4'), group.by = "SCT_snn_res.0.1") +
  ggtitle("CELL MARKER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Sub-clusters") +
  ylab("Features")
dev.copy(png, file = "240714-CLUSTER7-subset-dimplot-1059x2000.png", width = 1059, height = 2000)
dev.off()

# BY EXPRESSION DISTRIBUTION
FeaturePlot(cluster_7, features = c('TEK','PCDH12','SOX17','HEY1','DLL4',
                                    'RAMP2',
                                    'CD44','RUNX1','CD82','MYB','SPI1','NFE2','PTK2B','CORO1A'))
dev.copy(png, file = "240714-CLUSTER7-subset-expression-distribution-1000X1000.png", width = 1000, height = 1000)
dev.off()

DotPlot(cluster_7, features = c('TEK','PCDH12','SOX17','HEY1','DLL4',
                                'RAMP2',
                                'CD44','RUNX1','CD82','MYB','SPI1','NFE2','PTK2B','CORO1A'), group.by = "SCT_snn_res.0.1") +
  ggtitle("CELL MARKER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Sub-clusters") +
  ylab("Features")
dev.copy(png, file = "240714CLUSTER7-subset-expression-distribution-dotplot-800X800.png", width = 800, height = 800)
dev.off()



# ======================== # SUBSET CLUSTER 6 + 7 # ======================== # 
cluster6_7 <- readRDS(file = 'RDS/240714-cluster6+7(EHT+Blood)-subset.RDS')

cluster6_7 <- FindClusters(cluster6_7, dims= 1:15, resolution = 0.2)
c7 <- DimPlot(cluster6_7, reduction = "umap") 
original <- DimPlot(integrated, reduction = 'umap')
c7|original

de_markers6_7 <- FindAllMarkers(cluster6_7, only.pos = T) # for gene set ENRICHMENT analysis
filtered_markers <- de_markers6_7 %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, "240714_filtered(p)_markers_new-data-clusters6+7-subset.csv", row.names = FALSE)

filtered_markers <- read.csv(file = '240714_filtered(p)_markers_new-data-clusters6+7-subset-AUTOMATED.csv')

top10 <- filtered_markers  %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) # top 10/20 just for the heatmap 

# write.csv(top10, file ="240519-top10-avg_log2FC.csv")

DoHeatmap(cluster6_7, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240714-CLUSTER6+7-subset-DEGs-1059x705.png", width = 1059, height = 705)
dev.off()


# naive (immature, uncommitted progenitors): 'KIT', 'GATA2', + lack of lineage markers
# EHT: CD44','RUNX1','CD82'
# BLOOD LINEAGES
# Megakaryocyte:  'GP9', 'PF4', 'GATA1', 'TAL1', 'FLI1','ITGA2B','ICAM2','CD9'
# ERYTHROID: 'GYPA', 'KLF1', 'MYC', 'EPCAM','HBG1','HBG2' + 'HMGB1','HMGB2','TFAM','ZFP35L2','PA2G4','NOLC1','YBX1','PTTG1'
# GRANULOCYTE: 'AZU1', 'MPO','ITGB2','CEBP-D', 'CEBP-B', 'CEBP-A', 'CEBP-E', 
# myeloid: 'SPI1', 'MYB', 'GATA2'

# ENDOTHELIAL: 'TEK','PCDH12','SOX17','HEY1','DLL4',
# ENDOTHELIAL:'RAMP2',
# EHT: 'CD44','RUNX1','CD82',

# Myeloid (mega + ery + gran): PRTN3
# MEGA: ''TAL1', 'ITGA2B','ICAM2','CD9', (THEY'RE CLUSTER8)
# ERYTHROID: ''HBZ', 'GYPA','AHSP','HBE1' (THEY'RE CLUSTER9)
# GRANULOCYTE:'AZU1', 'MPO','ITGB2','CEBPD', 'CEBPB', 'CEBPA', 'RNASE2'
# LYMPHOID: 'HLA-DRA','HLA-DRB1')

FeaturePlot(cluster6_7, features = c('TEK','PCDH12','SOX17','HEY1','DLL4',
                                     'RAMP2',
                                     'CD44','RUNX1','CD82','MYB',
                                     'PRTN3',
                                     'TAL1', 
                                     'NFE2',
                                     'AZU1', 'MPO', 'RNASE2',
                                     'HLA-DRA','HLA-DRB1'))
dev.copy(png, file = '240714-Cluster6+7-cell-typing-3000x4000.png', width = 2000, heigh = 2000)
dev.off()



# ====
DotPlot(cluster6_7, features = c('TEK','PCDH12','SOX17','HEY1','DLL4',
                                 'RAMP2',
                                 'CD44','RUNX1','CD82','MYB',
                                 'PRTN3',
                                 'TAL1', 
                                 'NFE2',
                                 'AZU1', 'MPO', 'RNASE2',
                                 'HLA-DRA','HLA-DRB1'), group.by = "seurat_clusters") +
  ggtitle("CELL MARKER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Sub-clusters") +
  ylab("Features")
dev.copy(png, file = "240714-CLUSTER6+7-subset-expression-distribution-dotplot-800X800.png", width = 800, height = 800)
dev.off()

f1 <- DimPlot(cluster6_7, group.by = 'seurat_clusters', label = T, label.size = 6)
f2 <- DimPlot(integrated, group.by = 'SCT_snn_res.0.3', label = T, label.size = 6)

f1|f2
dev.copy(png, file = "2407014-cluster6+7-integrated-before-annotation-GO+cell-marker-UMAP-1500x708-new.png", width = 1500, height = 800)
dev.off()

# JUST FOR THE PPT
DimPlot(cluster6_7, split.by = 'seurat_clusters')
dev.copy(png, file = '240714-Cluster6+7-cell-typing-3000x4000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(integrated, features = c('TEK','PCDH12','SOX17','HEY1','DLL4',
                                     'RAMP2',
                                     'CD44','RUNX1','CD82','MYB',
                                     'PRTN3',
                                     'TAL1', 
                                     'NFE2',
                                     'AZU1', 'MPO', 'RNASE2',
                                     'HLA-DRA','HLA-DRB1'), split.by = 'Condition')
dev.copy(png, file = '240714-integrated-cell-typing-split-by-Condition-2000x2000.png', width = 2000, heigh = 6000)
dev.off()

#VASCULAR ENDOTHELIU,: 'NR2R2', 'KRT18','FLT4','SMARCA4','EPHB4', 'PCDH12'
#ARTERIAL ENDOTHELIUM: 'DLL4','GJA5','GJA4','EFNB2','SOX17','HEY1'
# ENDOTHELIUM: 'HOPX', 'TEK','
DotPlot(integrated, features = c('HOPX','TEK',
                                 'NR2F2','KRT18','PCDH12','FLT4','SMARCA4','EPHB4',
                                 'SOX17','HEY1','DLL4','GJA5','GJA4','EFNB2',
                                           'CD44','RUNX1','CD82','MYB',
                                           'PRTN3',
                                           'TAL1', 
                                           'NFE2',
                                           'AZU1', 'MPO', 'RNASE2',
                                           'HLA-DRA','HLA-DRB1'), group.by = 'Condition' )+
  ggtitle("GENE EXPRESSION ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Features") +
  ylab("Conditions")
dev.copy(png, file = "240714-intergrated-expression-distribution-dotplot-NEW2317-800X800.png", width = 800, height = 800)
dev.off()


FeaturePlot(integrated, features = c('TAL1', 'NFE2','AZU1', 'MPO', 'RNASE2','HLA-DRA','HLA-DRB1'), split.by = 'Condition')
dev.copy(png, file = '240714-integrated-cell-typing-blood-lineage-split-by-Condition-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(integrated, features = c('PRTN3', 'CD44','RUNX1','CD82','MYB'), split.by = 'Condition')
dev.copy(png, file = '240714-integrated-cell-typing-EHT-split-by-Condition-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(integrated, features = c('NR2F2','KRT18','PCDH12','FLT4','SMARCA4','EPHB4'), split.by = 'Condition')
dev.copy(png, file = '240714-integrated-cell-typing-vascular-endothelium-split-by-Condition-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(integrated, features = c('SOX17','HEY1','DLL4','GJA5','GJA4','EFNB2'), split.by = 'Condition')
dev.copy(png, file = '240714-integrated-cell-typing-arterial-endothelium-split-by-Condition-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

# ======================== # SUBSET CLUSTER 6, 8 # ======================== # 
de_markers_condition <- FindAllMarkers(cluster8_9, only.pos = T) # for gene set ENRICHMENT analysis
filtered_markers <- de_markers_condition %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, "240709_filtered(p)_markers_new-data-clusters8+9-subset.csv", row.names = FALSE)
#my.markers <- read.csv(file = '240709_filtered(p)_markers_new-data-clusters8+9-subset.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) # top 10/20 just for the heatmap 

# write.csv(top10, file ="240519-top10-avg_log2FC.csv")
sum(top10$gene %in% rownames(GetAssayData(cluster8_9, layer = 'scale.data')))
DoHeatmap(cluster8_9, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240709-CLUSTER8+9-subset-DEGs-1059x705.png", width = 1059, height = 705)
dev.off()

# ======================== # SUBSET CLUSTER 6, 7, 8, 9 # ======================== # 
de_markers_condition <- FindAllMarkers(cluster6789, only.pos = T)
filtered_markers_6789 <- de_markers_condition %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers_6789, "240709_filtered(p)_markers_new-data-clusters6789-subset.csv", row.names = FALSE)

top10_6789 <- filtered_markers_6789  %>% 
  group_by(cluster) %>% 
  top_n(20, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 

sum(top10_6789$gene %in% rownames(GetAssayData(cluster6789, layer = 'scale.data')))
DoHeatmap(cluster6789, features = top10_6789$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240709-CLUSTER6789-subset-DEGs-1059x705.png", width = 1059, height = 705)
dev.off()

# ======================== # SUBSET CLUSTER 6, 8 # ======================== # 
de_markers_condition <- FindAllMarkers(cluster6_8_blood, only.pos = T)
filtered_markers_6_8_blood <- de_markers_condition %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers_6_8_blood, "240709_filtered(p)_markers_new-data-clusters_6_8_blood-subset.csv", row.names = FALSE)

top10_6_8_blood<- filtered_markers_6_8_blood  %>% 
  group_by(cluster) %>% 
  top_n(20, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 

sum(top10_6_8_blood$gene %in% rownames(GetAssayData(cluster6_8_blood, layer = 'scale.data')))
DoHeatmap(cluster6_8_blood, features = top10_6_8_blood$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240709-cluster6_8_blood-subset-DEGs-1059x705.png", width = 1059, height = 705)
dev.off()

umap

# = try ANTO membrane marker for HSC-like cell =====
DotPlot(object = cluster6789, features = c("JAML", "CD44", "LCP1", "CYBA", 
                                           "AIF1", "CD34", "FABP5", "CORO1A", 
                                           "VAMP8", "ATP8B4", "PTPRC", "BST2", 
                                           "TNFSF10", "ANGPT1", "ITM2A", "VSIR", 
                                           "LY6E", "NKG7", "ITGA4", "NUCB2", "KATNB1",
                                           "CD53", "TAOK3", "TNFRSF1A", "CYTH4", "CD33",
                                           "ANKRD13D", "LAIR1", "LRPAP1", "MIF"), group.by = 'Condition' )+
  ggtitle("GENE EXPRESSION ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240709-test-dotplot-900x1500.png", width = 900, height = 1500)
dev.off() 

f1 <- FeaturePlot(integrated, features = c("JAML", "CD44", "LCP1", "CYBA", 
                                            "AIF1", "CD34", "FABP5", "CORO1A", 
                                            "VAMP8", "ATP8B4", "PTPRC", "BST2", 
                                            "TNFSF10", "ANGPT1", "ITM2A", "VSIR", 
                                            "LY6E", "NKG7", "ITGA4", "NUCB2", "KATNB1",
                                            "CD53", "TAOK3", "TNFRSF1A", "CYTH4", "CD33",
                                            "ANKRD13D", "LAIR1", "LRPAP1", "MIF"))
f1
dev.copy(png, file = "240709-test-integrated-dotplot-2000x4000.png", width = 2000, height = 4000)
dev.off() 
