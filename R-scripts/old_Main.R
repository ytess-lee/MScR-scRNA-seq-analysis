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
integrated_OG <- readRDS("RDS/210524-mechanobiology23_integretatedReps_050424.RDS") 
# no_prestim <- readRDS("RDS/240608-no-prestim.RDS")

# using new RDS file from Anto - 240614
#integrated <- readRDS("RDS/Mechanobiology_Harmonised_2reps_140624_TL.RDS") 

# FOR CLUSTER TOP 10
ElbowPlot(integrated_OG, ndims = 50)
# # # # # 
#  UMAP #
# # # # # 
integrated_OG <- FindNeighbors(integrated_OG, dims = 1:13)
integrated_OG <- FindClusters(integrated_OG, dims= 1:13, resolution = 0.4)
integrated_OG <- RunUMAP(object = integrated_OG, dims = 1:13) # it flips on the x-axis
# merged UAMP
DimPlot(integrated_OG, reduction = "umap")
dev.copy(png, file = "240705-CLUSTERS(res0.5)-OG-Data-937x708.png", width = 937, height = 708)
dev.off() 


# res = 0.4 =========================================================================================================================
deseq.res4 <- FindAllMarkers(integrated_OG,
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             only.pos = TRUE,
                             slot = 'counts')
write.csv(deseq.res4, file = '240705-findallmarker-across-clusters(res0.4)-OG-Data.csv')
filtered_markers <- deseq.res4 %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, file = '240705-findallmarker-across-clusters(res0.4)-OG-Data-filtered.csv')

top10 <- filtered_markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(top10, file = '240705-findallmarker-across-clusters(res0.4)-top10-filtered-OG-data.csv')

sum(top10$gene %in% rownames(GetAssayData(integrated_OG, layer = 'scale.data')))
DoHeatmap(integrated_OG, features = top10$gene) + theme(axis.text.y = element_text(size = 8))
dev.copy(png, file = "240705-CLUSTERS(res0.4)-top10-OG-Data-1300x900.png", width = 1300, height = 900)
dev.off()  # Close the PNG device






#saveRDS(no_prestim, "240608-no-prestim.RDS")
# # # # # # # # # # # # # # # # #
# DE GENES GROUPED BY CLUSTERS #
# # # # # # # # # # # # # # # # #
# temporarly store the cluster idents to be able to re-use them later
# object@slotName: to access a specific slot within an S4 object (OOP, not seurat specific)
temp_idents <- integrated_OG@active.ident

# assign the culture condition to the active identitites for the DEG analysis
Idents(integrated) <- "Condition" 

ElbowPlot(integrated, ndims = 50)

# # # # # # # # # # # 
# IMPORT AF MARKERS #
# # # # # # # # # # # 
all.markers <- read.csv(file = 'DEG_stimulation_merged_data(in).csv') # 3675
filtered.markers <- all.markers %>% filter(p_val_adj < 0.05) # 3655
# write.csv(filtered.markers, "240521_filtered_markers_AF.csv", row.names = FALSE)
filtered.markers <- read.csv(file ="240521_filtered_markers_AF.csv") # FROM ANTO

Top10 <- filtered.markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 

DoHeatmap(integrated, features = Top10$gene) + theme(axis.text.y = element_text(size = 6))



Top20 <- filtered.markers %>% 
  group_by(cluster) %>% 
  top_n(20, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 
# write.csv(Top20, "240526_Top20(AFcsv).csv", row.names = FALSE)
DoHeatmap(integrated, features = Top20$gene) + theme(axis.text.y = element_text(size = 6))

# # # # # # # # # # # # #
# CHECK GENE EXPRESSION #
# # # # # # # # # # # # #

# mechanical stimulus
FeaturePlot(object = integrated, features = c('TRPM1', 'TRPC6', 'PIEZO2')) # FROM Samidha
# osmotic stress
FeaturePlot(object = integrated, features = c('SLC4A11', 'KCNMA1', 'LRRC8A'))
# AQP <- related to vacuole formation during EHT
FeaturePlot(object = integrated, features = c('AQP1', 'AQP3', 'AQP5', 'AQP8', 'AQP9')) # low gene expression but doesn't mean low protein expression
FeaturePlot(object = integrated, features = c('PIEZO1', 'PIEZO2'))







# # # # # # # # # # # # ##
# FROM LITERAUTRE REVIEW #
# # # # # # # # # # # # ##
# CD326 <- Erythroid, NOT DETECTED
# CD9 <- Megakaryocyte, universal expression
# CD19 <- leukocyte, few dots at the HP lineage
FeaturePlot(object = integrated, features = c('CD9', 'CD19'))
FeaturePlot(object = integrated, features = c('YAP1', 'RUNX1', 'PG', 'HOXA', 'CD7', 'CD45'))


# https://cellregeneration.springeropen.com/articles/10.1186/s13619-024-00192-z
# GLUCOSE METABOLIC
# remove HK1 as it just express everywhere
# SLC2A1 & SLC2A3: highest expression in Laminar 
# HK1: less expression in pulsatile
FeaturePlot(object = integrated, features = c('SLC2A1', 'SLC2A3', 'HK2'))
VlnPlot(integrated, features = c('SLC2A1', 'SLC2A3', 'HK2'))

DotPlot(object = subset_seurat, features = c('SLC2A1', 'SLC2A3', 'HK2'), group.by = "Condition") +
  ggtitle("GLUCOSE METABOLIC") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()

# OXIDATIVE PHOSPHORYLATION - highly expressed across cluster
FeaturePlot(object = integrated, features = c('ATP5MC1', 'COX5A', 'CYC1'))
VlnPlot(integrated, features = c('ATP5MC1', 'COX5A', 'CYC1'))

DotPlot(object = subset_seurat, features = c('ATP5MC1', 'COX5A', 'CYC1'), group.by = "Condition") +
  ggtitle("OXIDATIVE PHOSPHORYLATION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()


# GLYCOLYSIS - highly expressed across cluster
FeaturePlot(object = integrated, features = c('LDHA', 'TPI1', 'PKM'))
VlnPlot(integrated, features = c('LDHA', 'TPI1', 'PKM'))

DotPlot(object = subset_seurat, features = c('LDHA', 'TPI1', 'PKM'), group.by = "Condition") +
  ggtitle("GLYCOLYSIS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()

# DETOXIFICATION OF ROS - # nothing too interesting
FeaturePlot(object = integrated, features = c('BNIP3', 'PRDX2', 'PRDX3', 'PRDX4', 'P4HB', 'GPX7'))
VlnPlot(integrated, features = c('BNIP3', 'PRDX2', 'PRDX3', 'PRDX4', 'P4HB', 'GPX7'))      
FeaturePlot(object = integrated, features = c('ATOX1', 'NUDT2', 'CYCS', 'TXN2'))
VlnPlot(integrated, features = c('ATOX1', 'NUDT2', 'CYCS', 'TXN2'))   


# FATTY ACID METABOLIC
# universally expressed across cluster, didn't save the img
FeaturePlot(object = integrated, features = c('SDHC', 'CA2', 'ADSL', 'ECHS1', 'DLD', 'ECI2'))
VlnPlot(integrated, features = c('SDHC', 'CA2', 'ADSL', 'ECHS1', 'DLD', 'ECI2'))      
# S100A10 exp (high to low): ENdo > lineage > bottleneck, saved img
FeaturePlot(object = integrated, features = c('ODC1', 'IDH1', 'S100A10', 'CCDC58'))
VlnPlot(integrated, features = c('ODC1', 'IDH1', 'S100A10', 'CCDC58'))  


# GLUTAMINE METABOLIC
# 'ASL', 'GOT1' GOT REMOVED DUE TO NON-SPECIFIC EXPRESSION
# FAH gpt upregulated in static and pulsatile
# PYCR1 specically up-regulated in all condition but prestim
FeaturePlot(object = integrated, features = c('GAD1', 'FAH', 'PYCR1', 'NIT2'))
VlnPlot(integrated, features = c('GAD1', 'FAH', 'PYCR1', 'NIT2'))      
# doesn't show much difference across cluster, img not saved
FeaturePlot(object = integrated, features = c('GOT2', 'ADSL', 'GPT2', 'GLUD1'))




# try Vln plot for DEG across clusters
VlnPlot(integrated, features = 'KLF2')
FeaturePlot(integrated, 'KLF2')

# try using other marker finding function but no luck
markers <- FindMarkers(object = integrated, ident.1 = 'Static', ident.2 = 'Pulsatile')
markers %>%
  dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) # uses the dplyr package to filter the results
write.csv(markers, "240529-Static-vs-Pulsatile.csv", row.names = TRUE)

FeaturePlot(integrated, features = c('ERC2-IT1', 'AC137770.1','LINC02398', 'AC097500.1',  'AC132153.1', 'AC104964.1'))
Top10 <- markers %>% 
  top_n(10, avg_log2FC) %>%
  arrange(desc(avg_log2FC)) 
write.csv(markers, "240529-Static-vs-Pulsatile.csv", row.names = TRUE)
svp.markers <- read.csv(file = '240529-Static-vs-Pulsatile.csv')
DoHeatmap(integrated, features = markers$gene) + theme(axis.text.y = element_text(size = 6))




# # # # # # # # #  
# PATHWAY GENE #
# # # # # # # # #  

# Reorder the levels of the condition factor
integrated$Condition <- factor(integrated_OG$Condition, levels = c("Prestimulation", "Static", "Pulsatile", "Laminar"))


# HIPPO/ YAP signalling
# tumor suppressor RASSF1A so it makes sense that it's not found
FeaturePlot(object = integrated, features = c('YAP1', 'TAZ', 'LATS1', 'LATS2', 'TEAD1'))

# VlnPlot to see which cluster has the highest expression (i doubt it)
VlnPlot(integrated, features = c('YAP1', 'TAZ', 'LATS1', 'LATS2', 'TEAD1'))

# Static/ laminar has higher TAZ and lower ANKRD1 expressions compared to other two
# I dont#t understand it :(
# = -0.03, no/ low correlation
FeatureScatter(integrated, 
               "TAZ", "ANKRD1", 
               group.by = "Condition", 
               pt.size = 1.5, 
               cols = condition.colour)

FeatureScatter(integrated, 
               "TAZ", "ANKRD1", 
               pt.size = 1.5, )

# Dot plot of gene expression across conditions
# really low expression on TAZ, removed to see if it scales better
DotPlot(object = integrated, features = c("YAP1", "LATS1", "LATS2", "TEAD1"), group.by = "Condition") +
  ggtitle("YAP Signalling Gene Expression in Different Conditions") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()

# let's try W/O prestim
# Subset to include only cells from the specified conditions
# subset_conditions <- c("Static", "Pulsatile", "Laminar")
# subset_seurat <- subset(integrated, subset = Condition %in% subset_conditions)

DotPlot(object = subset_seurat, features = c("YAP1", "LATS1", "LATS2", "TEAD1"), group.by = "Condition") +
  ggtitle("YAP Signalling Gene Expression in Selected Conditions") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()

# even though ShinyGO showed Ferroptosis in KEGG pathway analysis, 
# https://www.cell.com/cell/pdf/S0092-8674(23)00050-8.pdf
# Healthy human HSC: 'MYSM1', 'FTL', 'SLC40A1', 'SLC7A11', 'GPX4', 'ALAS1', 'RPL5'
# MYSM1 deficiency causes human HSC loss by ferroptosis, then we'll be expecting if pusaltile has an up-regulated ferroptosis, then expected low MYSM1 expression
# but MYSM1's got quite a spread on the featureplot, is there a way to find out what the value is caused it wasn't ion the csv anto sent
# consclusion: up-regulated ferroptosis but not MYSM1-deficiency dependent
DotPlot(object = integrated, features = c('MYSM1', 'FTL', 'SLC40A1', 'SLC7A11', 'GPX4', 'ALAS1', 'RPL5'), group.by = "Condition") +
  ggtitle("Ferroptosis-related Gene Expression in Different Conditions") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()

# 240605 update: these genes are not supposed to be highly expressed in ferroptosis-enriched HSPC/ HSCs
# unless the 'pathway enrichment' in ShinyGo is dependent on how many genes are up-regulated together, instead of just one? (Ferratin and SLC3A2)
# subset_conditions <- c("Static", "Pulsatile", "Laminar")
# subset_seurat <- subset(integrated, subset = Condition %in% subset_conditions)


DotPlot(object = subset_seurat, features = c('MYSM1', 'FTL', 'SLC40A1', 'SLC7A11', 'GPX4', 'ALAS1', 'RPL5', 'RPL24', 'RPS7', 'RPS18', 'NFE2L2'), group.by = "Condition") +
  ggtitle("Ferroptosis-related Gene Expression in Different Conditions") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()


FeaturePlot(object = subset_seurat, features = c('MYSM1', 'FTL', 'SLC40A1', 'SLC7A11', 'GPX4', 'ALAS1', 'RPL5', 'RPL24', 'RPS7', 'RPS18', 'NFE2L2'))



# https://www.nature.com/articles/s41422-020-00441-1/figures/5
# GHS production '
FeaturePlot(object = subset_seurat, features = c('SLC7A11', 'GSR', 'GSS', 'GCLM', 'GCLC'))

# ROS DETOXIFICATION
FeaturePlot(object = subset_seurat, features = c('GPX4', 'NQO1', 'HMOX1', 'ALDH3A1', 'AKR1C'))

# IRON METABOLISM
FeaturePlot(object = subset_seurat, features = c('MT1G', 'FTL', 'FTH1', 'SLC40A1', 'FECH'))




# # # # # # # # # # 
# Co-localisation #
# # # # # # # # # # 

# Define a named vector of colours, same colours as thowed on the heatmap
condition.colour <- c('Prestimulation' = '#F8766D',
                      'Static' = '#C2D989',
                      'Laminar' = '#C476FF',
                      'Pulsatile' = '#00BFC4')

# Can't change how the conditions are arranged in the figure legend, ggplot maybe?
# co-localisation by conditions
FeatureScatter(integrated, 
               "PIEZO1", "PIEZO2", 
               group.by = "Condition", 
               pt.size = 1.5, 
               cols = condition.colour)

# co-localisation in one pot
FeatureScatter(integrated, 
               "PIEZO1", "PIEZO2",
               pt.size = 1.5,
               cols = condition.colour) +
  facet_wrap(Idents(integrated))

















# export gene list for the heatmap generated
# write.csv(Top10, "240521_Top10(AFcsv).csv", row.names = FALSE)

# # # # # # # #
# CELL TYPING #
# # # # # # # #
# ENDOTHELIAL CELL MARKER
# https://www.proteinatlas.org/ENSG00000187513-GJA4/tissue+cell+type
FeaturePlot(object = integrated, features = c('DLL4','GJA5', 'EFNB2', 'HEY1', 'NOTCH4', 'SOX17', 'EPHB4', 'NR2F2', 'FLT4'))
# ARTERIAL: 'DLL4', 'GJA4', 'GJA5', 'EFNB2', 'SOX17', 'HEY1'
# HEY2 is more smooth muscle related
# GJA4 leaked passed bottleneck
# FeaturePlot(object = integrated, features = c('DLL4', 'GJA4', 'GJA5', 'EFNB2', 'SOX17', 'HEY1'))
FeaturePlot(object = integrated, features = c('DLL4','GJA4', 'GJA5', 'EFNB2', 'HEY1', 'NOTCH4', 'SOX17'))

# Venous: 'NR2F2', 'EPHB4', 'SMARCA4', 'FLT4'
# a bit leakage of EMCN towards the bottleneck and PASSED the bottleneck so removed
# FeaturePlot(object = integrated, features = c('NR2F2', 'EPHB4', 'FLT4'))
FeaturePlot(object = integrated, features = c('EPHB4', 'NR2F2', 'FLT4'))

# EHT marker: 'RUNX1', 'SOX17', 'GATA2', 'MEIS1', 'GFI1', 'CD44' (CD41, CD43, CD45,  'CD31', 'VEC', 'FLK-1' not founc)
# GFI1/ GFI1B: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9136496/
# MEIS1 (earliest EHT marker + prior to RUNX1): https://pubmed.ncbi.nlm.nih.gov/37500618/
# EHT markers: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6490701/
# I don't understand what's causing no/low SOX17 expression pass the endothelial cells
FeaturePlot(object = integrated, features = c('RUNX1', 'SOX17', 'GATA2', 'MEIS1', 'GFI1', 'GFI1B', 'CD41', 'CD31', 'VEC', 'FLK-1'))


# Haematopoietic lineage markers: 'KLF1', 'CD14', 'ITGA2B', 'MPO', 'CSF2RA', 'CD14' 
# 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'CD1A' 
# klf1, ITGA2B: erythoid
# righty
FeaturePlot(object = integrated, features = c('KLF1', 'ITGA2B'))


# lymphoid: https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.b.21614
# T, B, NK CELLS <- CD3, CD19, CD56, HLA-DR, CD34 
# CD3, CD56 NOT FOUND
# cd19: super low so removed
# right-righty
FeaturePlot(object = integrated, features = c('HLA-DRA', 'HLA-DRB1', 'HLA-DRB5'))


# myeloid: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4742930/; https://www.biocompare.com/Editorial-Articles/612866-A-Guide-to-Myeloid-Cell-Markers/
# not found: CD45, CD206, CD169, CD11c, CD123, CD16, CD15, Ly6C, CD1a, CD141, CD11b, CD64, CD1c, CD103, CD11b
# NON-SPECIFIC PASS BOTTLENECK ARE EXCLUDED: 'ANGPT2','ARG1', 'CCL2','CSF1','ANPEP', 'ADGRE1','CD24', IL3RA
# LOW EXPRESSION: CD1A
# IL3RA/ CD123: responsible for the proliferation, survival, and differentiation of hematopoietic cells.
# MPO: peroxidase enzyme
# CSF2RA: granulocytes and macrophages
# righy-lefty + MPO lol
FeaturePlot(object = integrated, features = c( 'CX3CR1', 'CSF2RA', 'CD14', 'CCR2', 'CD1A', 'MPO'))

# HAEMATOPOEITIC LINEAGE MASH
FeaturePlot(object = integrated, features = c( 'CX3CR1', 'CSF2RA', 'CD14', 'CCR2', 'CD1A', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'KLF1', 'ITGA2B', 'MPO'))








# ============================================================================================== #

# # # # # # # # #
# USED SCRIPTS  #
# # # # # # # # #

# for subsetting
data1 <- subset(integrated, subset = Condition == c("Prestimulation","Static"))


# determine the # of dimensions/ PCs we need to capture the majority of the variation in the integrated
ElbowPlot(integrated_OG, ndims = 50)
# # # # # 
#  UMAP #
# # # # # 
integrated_OG <- FindNeighbors(integrated_OG, dims = 1:13)
integrated_OG <- FindClusters(integrated_OG, dims= 1:13, resolution = 0.5)
integrated_OG <- RunUMAP(object = integrated_OG, dims = 1:13) # it flips on the x-axis
# merged UAMP
DimPlot(integrated_OG, reduction = "umap")
# split UMAP by conditions
DimPlot(integrated_OG, reduction = "umap",
        split.by = "dataset", group.by = "Condition",ncol = 2) # number of columns for display when combining plots


# # # # # 
# tSNE  # - 240606 Didn't save tSNE coz it doesn't say much about the clusters and it's ugly lol
# # # # #
data <- FindNeighbors(object = integrated, dims = 1:13, verbose = FALSE)
data <- RunTSNE(data, dims = 1:13)
# merged tSNE
DimPlot(data, reduction ="tsne")
# split tSNE by conditions
DimPlot(data, reduction = "tsne",
        split.by = "dataset", group.by = "Condition", ncol = 2)

# # # # # #
# TABLES  #
# # # # # #
table(integrated@active.ident)
table(Idents(integrated))

# # # # # # # # # # # # # # # # # # # # ## # # # # 
# FILTERED MARKERS W/ p_val_adj < 0.05 && EXPORT #
# # # # # # # # # # # # # # # # # # # # ## # # # # 
# Apply SCTransform normalization
integrated <- SCTransform(integrated, verbose = TRUE)

# Prepare SCT object for marker detection
integrated <- PrepSCTFindMarkers(integrated)

de_markers_condition <- FindAllMarkers(integrated_OG, only.pos = T) # for gene set ENRICHMENT analysis
filtered_markers <- de_markers_condition %>% filter(p_val_adj < 0.05)
write.csv(filtered_markers, "240629_filtered(p)_markers_OG-DATA.csv", row.names = FALSE)
my.markers <- read.csv(file = '240629_filtered(p)_markers_OG-DATA.csv')

# # # # # # # # # # # # # # # # # # #
# TOP10 (10 genes for EACH CLUSTER) #
# # # # # # # # # # # # # # # # # # #
top10_og_func <- my.markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) # top 10/20 just for the heatmap 

# write.csv(top10, file ="240519-top10-avg_log2FC.csv")
sum(top10_og_func$gene %in% rownames(GetAssayData(integrated_OG, layer = 'scale.data')))
DoHeatmap(integrated_OG, features = top10_og_func$gene) + theme(axis.text.y = element_text(size = 6))

# MINE MARKER
# Using the dataset recieved from the memory stick, named integrated
# filtered out p>0.05 and use top10 to select the most sign based on avg_log2FC
# exported the csv and do heat map using that csv
# 15 genes got omitted as they were not found in the scale.data slot for the SCT assay
my.markers <- read.csv(file = '240521_Marker_Clusters_merged.csv')
my.filtered.markers <- my.markers %>% filter(p_val_adj < 0.05)


my.top10 <- my.filtered.markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) 
# write.csv(my.top10, "240526_Top10(OG, filtered).csv", row.names = FALSE)

DoHeatmap(integrated, features = my.top10$gene) + theme(axis.text.y = element_text(size = 6))
FeaturePlot(object = integrated, features = c('FAM167B', 'A2M', 'KLF2', 'CCL4L2'))



# # # # # # # # # # # # # # # # # # # # ## # # # # # # # # #
# # Check gene marker expression distribution in clusters #
# # # # # # # # # # # # # # # # # # # # ## # # # # # # # # 
# # Prestimulation: "ADAMTS18", "SOST", "REN", "CCN4", "DCN", "ADAMTS20", "AC091178.2", "BMP10", "WNT2", "SST"
FeaturePlot(object = integrated, features = c("ADAMTS18", "SOST", "REN", "CCN4", "DCN", "ADAMTS20"))
FeaturePlot(object = integrated, features = c("AC091178.2", "BMP10", "WNT2", "SST"))

# # static: NA

# # pulsatile: CCL4L2, CRYAB, POU3F1, FNDC7
FeaturePlot(object = integrated, features = c("CCL4L2", "CRYAB", "POU3F1", "FNDC7", "AC020651.2", "LINC00880"))

# # laminar: KLF2, COL15A1, NCALD, A2M, FAM167B, IGFBP5
FeaturePlot(object = integrated, features = c("KLF2", "COL15A1", "NCALD", "A2M", "FAM167B", "IGFBP5"))


# # # # # # # # # # # # # # # # # # # # #
# CHECK HEATMAP WARNING MSG EXPRESSION  #
# # # # # # # # # # # # # # # # # # # # #

# filtered top 10 && ranked by avg_log2FC (that were omitted)
#1 'Z99289.2', 'ADAMTSL2', 'LYPD5', 'Z99289.1', 'HSPA12B', 'AC108134.1' 
FeaturePlot(object = integrated, features = c('Z99289.2', 'ADAMTSL2', 'LYPD5', 'Z99289.1', 'HSPA12B', 'AC108134.1'))
#2 'S100A14', 'AMPD1', 'AC100771.2', 'AC078860.2', 'NKAPL' 
FeaturePlot(object = integrated, features = c('S100A14', 'AMPD1', 'AC100771.2', 'AC078860.2', 'NKAPL'))
#2 'TTR', 'Z93241.1', 'AC239799.2', 'AC012447.1'
FeaturePlot(object = integrated, features = c('TTR', 'Z93241.1', 'AC239799.2', 'AC012447.1'))




# filtered top 10 && ranked by p_val_adj (that were omitted)
#1 'AC002519.1', 'LRRIQ3', 'AL080317.1', 'LIME1', 'CIART', 'LINC00211'
FeaturePlot(object = integrated, features = c('AC002519.1', 'LRRIQ3', 'AL080317.1', 'LIME1', 'CIART', 'LINC00211'))
#2 'TRIM34', 'NMNAT2', 'OOEP', 'SOS1-IT1', 'MVP', 'SLFN5', 
FeaturePlot(object = integrated, features = c('TRIM34', 'NMNAT2', 'OOEP', 'SOS1-IT1', 'MVP', 'SLFN5'))
#3'AL139021.2', 'TJP3', 'PP2D1', 'SLC4A3', 'LINC01128', 'CPEB2-DT', '
FeaturePlot(object = integrated, features = c('AL139021.2', 'TJP3', 'PP2D1', 'SLC4A3', 'LINC01128', 'CPEB2-DT'))
#4'SLC5A10', 'AC010168.2', 'AC022973.4', 'CASC9', 'ZNF461', 'AP000845.1'
FeaturePlot(object = integrated, features = c('SLC5A10', 'AC010168.2', 'AC022973.4', 'CASC9', 'ZNF461', 'AP000845.1'))
#5'DTX2', 'PDE8A', 'LINC02742', 'NLGN2', 'AP000766.1', 'AC005280.2'
FeaturePlot(object = integrated, features = c('DTX2', 'PDE8A', 'LINC02742', 'NLGN2', 'AP000766.1', 'AC005280.2'))
#6'NT5DC3', 'AC093298.2', 'CYB561', 'ANKRD40', 'AC007773.1', 'MCM3AP', 
FeaturePlot(object = integrated, features = c('NT5DC3', 'AC093298.2', 'CYB561', 'ANKRD40', 'AC007773.1', 'MCM3AP'))
#7 'ENPP5', 'MMGT1'
FeaturePlot(object = integrated, features = c('ENPP5', 'MMGT1'))


sum(top10$gene %in% rownames(GetAssayData(data, layer = 'scale.data')))
[1] 25

sum(filtered.markers$gene %in% rownames(GetAssayData(data, layer = 'scale.data')))
[1] 3510
# 3510 of them have corresponding entries in the scaled data matrix of the Seurat object

# # # # # # # # # # # 
# WARNING MESSAGES  #
# # # # # # # # # # # 
# P_VAL_ADJ
# PRESTIM
FeaturePlot(object = integrated, features = c('CYB561', 'MCM3AP'))
# STATIC
FeaturePlot(object = integrated, features = c('MVP', 'CASC9', 'NLGN2', 'DTX2', 'ZNF461'))
# LAMINAR
FeaturePlot(object = integrated, features = c('LINC01128', 'AL139021.2', 'PDE8A', 'ANKRD40', 'MMGT1'))
# PULSATILE
FeaturePlot(object = integrated, features = c('SLC4A3', 'SLFN5'))

# # # # # # # # # # # # # 
# CHECK GENE EXPRESSION # After top10 DE analysis using AF-csv
# # # # # # # # # # # # #
# PRESTIM
# 'GAL', 'CKS1B', 'EDN1', 'KRT18', 'DIAPH3',	'GPC6',	'KRT8',	'CENPK',	'CEP128',	'CENPW'
FeaturePlot(object = integrated, features = c('GAL',	'CKS1B',	'EDN1',	'KRT18'))
FeaturePlot(object = integrated, features = c('DIAPH3',	'GPC6',	'KRT8',	'CENPK'))
FeaturePlot(object = integrated, features = c('CEP128',	'CENPW'))

# STATIC
# 'RGS6'	'ACSM3'	'PRTN3'	'PLEK'	'LYZ'	'SRGN'	'LTBP1'	'PLXDC2'	'S100A8'	'SLC24A3'
FeaturePlot(object = integrated, features = c('RGS6', 'ACSM3', 'PRTN3', 'PLEK'))
FeaturePlot(object = integrated, features = c('LYZ',	'SRGN',	'LTBP1',	'PLXDC2'))
FeaturePlot(object = integrated, features = c('S100A8', 'SLC24A3'))

# PULSATILE 240521-top10-PULSATILE-AFcsv-1
#'HBG2',	'HBG1',	'HBA1',	'HBZ'	'HBA2',	'S100A4',	'PPBP',	'HBD',	'S100B',	'MT2A'
FeaturePlot(object = integrated, features = c('HBG2',	'HBG1',	'HBA1',	'HBZ'))
FeaturePlot(object = integrated, features = c('HBA2',	'S100A4',	'PPBP',	'HBD'))
FeaturePlot(object = integrated, features = c('S100B',	'MT2A'))


# LAMINAR 240521-top10-LAMINAR-AFcsv-1
#'IGFBP5'	'NQO1'	'LAMA4'	'SLC9A3R2'	'COL15A1'	'A2M'	'HECW2'	'PLVAP'	'CD34'	'LBH'
FeaturePlot(object = integrated, features = c('IGFBP5',	'NQO1',	'LAMA4',	'SLC9A3R2'))
FeaturePlot(object = integrated, features = c('COL15A1',	'A2M',	'HECW2',	'PLVAP'))	
FeaturePlot(object = integrated, features = c('CD34',	'LBH'))
