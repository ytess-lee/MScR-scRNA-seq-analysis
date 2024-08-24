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
# remove: PECA<1, KDR, COL4A2
FeaturePlot(integrated, features = c('TEK','PCDH12', 'RAMP2', 'HEY1','HOPX', 'MECOM','DLL4', 'KRT18', 'NR2F2', 'SOX17'))
dev.copy(png, file = "240707-Endothelial-900X700.png", width = 900, height = 700)
dev.off()


# Calcium binding activity
FeaturePlot(integrated, features = c('MCTP1', 'MCTP2', 'TRPV4'))




# dotplot to visualise the distribution of gene experssion across clusters
DotPlot(object = integrated, features = c('NFE2', 'DOK2', 'PTK2B', 'CORO1A', 'SPI1' , 'MYB', 'KLF1', 'GATA1', 'GFI1',
                                          'GATA2', 'ITGA2B', 'CD44' ,'RUNX1' ,'CD82',
                                          'TEK','PCDH12', 'RAMP2', 'HEY1','HOPX', 'MECOM','DLL4', 'KRT18', 'NR2F2'), group.by = "Condition") +
  ggtitle("GENE EXPRESSION ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")
dev.copy(png, file = "240707-Endothelial+EHT+lineages-dotplot-900x705.png", width = 900, height = 705)
dev.off() 

# = = = = = = = = = = = = = = = = = AF 20220 BLOOD MARKERS = = = = = = = = = = = = = = = = = = = = = = = = = = #

# Mono-DC
# "IGLL1", "GZMA", "KLRB1", "DNASE1L3", "IGKC", "C15orf48", "LYVE1", "FCER1A", "CDA", "CCL5", "AZU1", "LILRB5", "LGALS2", "CXCL3", "CXCL2", "CETP", "AGR2", "NR4A3", "TIMD4", "IFITM1"
# "OSM", "ABCA1", "S100A8", "LYZ", "HLA-DRA", "S100A12", "CKAP4", "VMO1", "G0S2", "S100A9", "CD163", "HLA-DPA1", "MPO", "S100P", "MS4A3", "HLA-DPB1", "CXCL8", "ATP6V0C", "CCL4L2", "CXCL16" 
# "FOLR2", "RNASE1", "CD74", "HLA-DRB1", "LYST", "PADI4", "RETN", "C1orf54", "C5AR1", "RNASE2", "TSPO", "CLEC5A", "AHSP", "SNAI1", "HBM", "MARCO", "AREG", "SPIB", "PLA2G7", 
# "MCEMP1", "RGS10", "AVPI1", "CLEC9A", "HLA-DMA", "CST3", "MAF", "FCN1", "MT1G", "CLEC11A", "CFD", "HLA-DMB", "BAG3", "AXL", "MGST1", "HP", "ZFP36", "ANKRD37", "MALAT1", "SERPINB1" 
# "IL1B", "TNFAIP3", "PTGER4", "ALB", "APLP2", "STXBP2", "CSTA", "EPX", "MSLN", "FOS", "HSPE1", "RASGRP2", "DNAJB4", "HCST", "MCOLN1", "HSPA6", "DDIT4", "PRSS57", "MEF2C", "FTH1"
# "FGL2", "MT-ATP8", "CCL4", "CD7", "IFI16", "HBEGF", "MT-CO1", "MSRB1", "PTPRC", "HBG1", "ID2", "RGS2", "JUNB", "B2M", "HIST1H1E", "GATM", "GLIPR1", "IGFBP2", "ATP2B1", "PSMB9"
# "SAT1", "NME2", "ITM2B", "JAML", "CCL3", "ARRDC3", "CTSL", "LGALS1", "LTB", "HLA-DQB1", "S100A10", "CTSS", "MT-ND5", "APOE", "RGS1", "SGK1", "HLA-DRB5", "IFITM2", "HERPUD1", "SPINT2"
# "HSPD1", "ZFAND2A", "NEAT1", "MNDA", "ANXA2", "MS4A7", "DUSP2", "GLUL", "EVL", "DUSP1", "ABI3", "EZR", "PLIN2", "S100B", "HAVCR2", "ALDH2", "TUBA1A", "UBC", "RHOB", "LPAR6"
# "CTSB", "HSPA8","TIMP1", "ARL4C", "CRIP1", "LAPTM4A", "NFKBIZ", "NR4A2", "CYTIP", "HSPH1", "HSPA1A", "SQSTM1", "CD4", "RNASE3", "HSPA1B", "AHNAK", "FTL", "KCNK6", "PPP1R15A", "ANXA1"
# "HLA-DQA1", "MCL1", "SH3BP5", "NR4A1", "RNASE6", "PPP1R14A", "KLF6", "REL", "CD86"

# TO INCLUDE IN DISSO
# "FCER1A", "AGR", 'HLA-DRA', "HLA-DPA1", "HLA-DPB1", "CD74", "HLA-DRB1", "HLA-DMA",
# "CDA", "AZU1", "S100A9", "MPO","S100P", "MS4A3", "LYST", "RETN", "MCEMP1", "FCN1",
# "PLA2G7",
# "S100AB", 'LYZ',"CFD", "FGL2"

# 20
FeaturePlot(object = integrated, features = c("IGLL1", "GZMA", "KLRB1", "DNASE1L3", "IGKC", "C15orf48", "LYVE1", "FCER1A", "CDA", "CCL5", "AZU1", "LILRB5", "LGALS2", "CXCL3", "CXCL2", "CETP", "AGR2", "NR4A3", "TIMD4", "IFITM1"))
dev.copy(png, file = "240709-Mono-DC-10-20-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 40
FeaturePlot(object = integrated, features = c("OSM", "ABCA1", "S100A8", "S100A12", "CKAP4", "VMO1", "G0S2", "S100A9", "CD163", "HLA-DPA1", "MPO", "S100P", "MS4A3", "HLA-DPB1", "CXCL8", "ATP6V0C", "CCL4L2", "CXCL16"))
dev.copy(png, file = "240709-Mono-DC-20-40-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 60
FeaturePlot(object = integrated, features = c("FOLR2", "RNASE1", "CD74", "HLA-DRB1", "LYST", "PADI4", "RETN", "C1orf54", "C5AR1", "RNASE2", "TSPO", "CLEC5A", "AHSP", "SNAI1", "HBM", "MARCO", "AREG", "SPIB", "PLA2G7"))
dev.copy(png, file = "240709-Mono-DC-40-60-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 80
FeaturePlot(object = integrated, features = c("MCEMP1", "RGS10", "AVPI1", "CLEC9A", "HLA-DMA", "CST3", "MAF", "FCN1", "MT1G", "CLEC11A", "CFD", "HLA-DMB", "BAG3", "AXL", "MGST1", "HP", "ZFP36", "ANKRD37", "MALAT1", "SERPINB1"))
dev.copy(png, file = "240709-Mono-DC-60-80-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 100
FeaturePlot(object = integrated, features = c("IL1B", "TNFAIP3", "PTGER4", "ALB", "APLP2", "STXBP2", "CSTA", "EPX", "MSLN", "FOS", "HSPE1", "RASGRP2", "DNAJB4", "HCST", "MCOLN1", "HSPA6", "DDIT4", "PRSS57", "MEF2C", "FTH1"))
dev.copy(png, file = "240709-Mono-DC-80-100-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 120
FeaturePlot(object = integrated, features = c("FGL2", "MT-ATP8", "CCL4", "CD7", "IFI16", "HBEGF", "MT-CO1", "MSRB1", "PTPRC", "HBG1", "ID2", "RGS2", "JUNB", "B2M", "HIST1H1E", "GATM", "GLIPR1", "IGFBP2", "ATP2B1", "PSMB9"))
dev.copy(png, file = "240709-Mono-DC-100-120-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 140
FeaturePlot(object = integrated, features = c("SAT1", "NME2", "ITM2B", "JAML", "CCL3", "ARRDC3", "CTSL", "LGALS1", "LTB", "HLA-DQB1", "S100A10", "CTSS", "MT-ND5", "APOE", "RGS1", "SGK1", "HLA-DRB5", "IFITM2", "HERPUD1", "SPINT2"))
dev.copy(png, file = "240709-Mono-DC-120-140-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 160
FeaturePlot(object = integrated, features = c("HSPD1", "ZFAND2A", "NEAT1", "MNDA", "ANXA2", "MS4A7", "DUSP2", "GLUL", "EVL", "DUSP1", "ABI3", "EZR", "PLIN2", "S100B", "HAVCR2", "ALDH2", "TUBA1A", "UBC", "RHOB", "LPAR6"))
dev.copy(png, file = "240709-Mono-DC-140-160-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 180
FeaturePlot(object = integrated, features = c("CTSB", "HSPA8","TIMP1", "ARL4C", "CRIP1", "LAPTM4A", "NFKBIZ", "NR4A2", "CYTIP", "HSPH1", "HSPA1A", "SQSTM1", "CD4", "RNASE3", "HSPA1B", "AHNAK", "FTL", "KCNK6", "PPP1R15A", "ANXA1", "HLA-DQA1", "MCL1", "SH3BP5", "NR4A1", "RNASE6", "PPP1R14A", "KLF6", "REL", "CD86"))
dev.copy(png, file = "240709-Mono-DC-180-xxx-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# THE ONE TO INCLUDE FOR CELL TYPING
FeaturePlot(object = integrated, features = c("FGL2", 'HLA-DRA', "HLA-DPA1", "HLA-DPB1", "CD74", "HLA-DRB1", "HLA-DMA",
                                              "AZU1", "S100A8", "MPO","MS4A3", "LYST", "RETN",  "FCN1",
                                              "PLA2G7", "CFD","S100A9", "LYZ"))
dev.copy(png, file = "240709-Mono-DC-2000x2000.png", width = 2000, height = 2000)
dev.off() 


# Neutrophil Myeloid Progenitors
# "PRG2", "CCL4", "CLC", "CCL3", "CXCL8", "CST7", "MNDA", "CSTA", "HSPA6", "ALDH1A1", "CXCL2", "RETN", "CTSG", "HLA-DRB5", "LYZ", "HBA2", "AZU1", "NFKBIZ", "HBG2", "MPO", "PRTN3", "MS4A3", "HBA1", "RAMP1", "ALB", "SRGN", "TSPO", "RNASE2", "S100A9", "SSR4", "PRDX1", "CEBPD", "APOA2", "S100A11", "MZB1", "RNASE3", "MALAT1", "BTG2", "APOA1", "SPINK2"
# "NKG7", "ZFAND2A", "CD37", "FTH1", "PNP", "TMSB4X", "SLC4A1", "SDF2L1", "CFD", "CD52", "HSP90B1", "HBM", "PRDX2", "MEST", "AHSP", "MDK", "CLEC11A", "GATA2", "AFP", "SOX4", "FLNA", "SERPINA1", "APLP2", "S100P", "MT1G", "SPINT2", "HMGA1", "SMIM24", "RGS10", "AC104389.1", "GYPA", "CD99", "CD74", "ICAM3", "STXBP2", "LSP1", "PRSS57", "FKBP4", "LY6E", "MEG3"
# "TUBB", "HCST", "IL1B", "PLAC8", "S100A8", "TYROBP", "FKBP1A", "SERPINB1", "CYTL1", "CELF2", "IFITM2", "HSPD1", "TRGC1", "HBG1", "MARCKSL1", "HOXA9", "EGFL7", "HEMGN", "TUBA1A", "DDAH2", "HLA-DRA", "TAOK3", "BST2", "CD34", "HPGDS", "MT1H", "ISG15", "NRIP1", "ALAS2", "HLA-DMA", "LAPTM4A", "HBZ", "HBB", "UROD", "PMAIP1", "DEFA4"

# 40
FeaturePlot(object = integrated, features = c("PRG2", "CLC", "CCL3", "CXCL8", "CST7", "MNDA", "CSTA", "HSPA6", "ALDH1A1", "CXCL2", "RETN", "CTSG", "HLA-DRB5", "LYZ", "HBA2", "AZU1", "NFKBIZ", "HBG2", "MPO", "PRTN3", "MS4A3", "HBA1", "RAMP1", "ALB", "SRGN", "TSPO", "RNASE2", "S100A9", "SSR4", "PRDX1", "CEBPD", "APOA2", "S100A11", "MZB1", "RNASE3", "MALAT1", "BTG2", "APOA1", "SPINK2"))
dev.copy(png, file = "240709-Neutrophil Myeloid Progenitors-40-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 80
FeaturePlot(object = integrated, features = c("NKG7", "ZFAND2A", "CD37", "FTH1", "PNP", "TMSB4X", "SLC4A1", "SDF2L1", "CFD", "CD52", "HSP90B1", "HBM", "PRDX2", "MEST", "AHSP", "MDK", "CLEC11A", "GATA2", "AFP", "SOX4", "FLNA", "SERPINA1", "APLP2", "S100P", "MT1G", "SPINT2", "HMGA1", "SMIM24", "RGS10", "AC104389.1", "GYPA", "CD99", "CD74", "ICAM3", "STXBP2", "LSP1", "PRSS57", "FKBP4", "LY6E", "MEG3"))
dev.copy(png, file = "240709-Neutrophil Myeloid Progenitors-80-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 120
FeaturePlot(object = integrated, features = c("TUBB", "HCST", "IL1B", "PLAC8", "S100A8", "TYROBP", "FKBP1A", "SERPINB1", "CYTL1", "CELF2", "IFITM2", "HSPD1", "TRGC1", "HBG1", "MARCKSL1", "HOXA9", "EGFL7", "HEMGN", "TUBA1A", "DDAH2", "HLA-DRA", "TAOK3", "BST2", "CD34", "HPGDS", "MT1H", "ISG15", "NRIP1", "ALAS2", "HLA-DMA", "LAPTM4A", "HBZ", "HBB", "UROD", "PMAIP1", "DEFA4"))
dev.copy(png, file = "240709-Neutrophil Myeloid Progenitors-120-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# 10: MNDA, CSTA, PRTN3, SRGN, CEBPD, NKG7, CD52, LSP1, PLAC8, TYROBP
FeaturePlot(object = integrated, features = c('MNDA', 'CSTA', 'PRTN3', 'SRGN', 'CEBPD', 'NKG7','CD52', 'LSP1', 'PLAC8', 'TYROBP',
                                              "TUBB", "HCST", "IL1B", "PLAC8", "S100A8", "TYROBP", "FKBP1A", "SERPINB1"))
dev.copy(png, file = "240709-Neutrophil Myeloid Progenitors-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# MONOCYTE
# "XCL1", "EREG", "GZMA", "IGLL1", "AREG", "CCL20", "KLRB1", "IFI27", "C15orf48", "IL32", "CD7", "VCAM1", "FABP3", "RARRES3", "CLEC2B", "CXCL12", "RND3", "A2M", "MSR1", "ISG20", "BIRC3", "ZNF331",  "ITGA4", "C2", "SDC3",   "CLEC10A",  "MAN2B1", "PLBD1", "S100A6", "SELL", "RNF144B", "PTGS2",   "GPR183",  "PRDM1",  "HES1",  "RBP7", "CAPG"
# "MSRB1", "MS4A4A", "CD69", "S100A4", "VIM", "TNF", "VCAN", "LTA4H", "LPCAT1","EPB41L2", "NCF1", "EVI2B", "PFKFB3", "ATP1B1", "CMC1", "HSP90AA1", "HMOX1", "C1QA",  "LIPA", "JAK1" , "LDHA", "FCGRT", "ICAM1", "HSP90AB1", "TSPAN4", "VAMP5", "CHI3L1", "MT-CO2", "MT-ND4L", "GADD45B", "MIF", "IFITM3", "CYP1B1"
# "THBS1", "PHACTR1", "PYCARD", "MS4A7", "MAP3K8", "DUSP2", "ID1", "MTSS1", "SMPDL3A", "ATF3",   "PTMA", "NPC2", "CCL2", "CTSC",  "EGR1", "IFNGR1", "KLF2", "CPVL",  "MRC1",   "TCF4", "UTRN", "CDKN1A", "PDE4B", "CALR", "ZFP36L1", "PDK4", "IFIT2",  "TPP1",  "BHLHE40", "AZU1",   "BNIP3L", "ID3"

FeaturePlot(object = integrated, features = c("XCL1", "EREG", "GZMA", "IGLL1", "AREG", "CCL20", "KLRB1", "IFI27", "C15orf48", "IL32", "CD7", "VCAM1", "FABP3", "RARRES3", "CLEC2B", "CXCL12", "RND3", "A2M", "MSR1", "ISG20", "BIRC3", "ZNF331",  "ITGA4", "C2", "SDC3",   "CLEC10A",  "MAN2B1", "PLBD1", "S100A6", "SELL", "RNF144B", "PTGS2",   "GPR183",  "PRDM1",  "HES1",  "RBP7", "CAPG"))
dev.copy(png, file = "240709-Mono-80-2000x2000.png", width = 2000, height = 2000)
dev.off() 

FeaturePlot(object = integrated, features = c("MSRB1", "MS4A4A", "CD69", "S100A4", "VIM", "TNF", "VCAN", "LTA4H", "LPCAT1","EPB41L2", "NCF1", "EVI2B", "PFKFB3", "ATP1B1", "CMC1", "HSP90AA1", "HMOX1", "C1QA",  "LIPA", "JAK1" , "LDHA", "FCGRT", "ICAM1", "HSP90AB1", "TSPAN4", "VAMP5", "CHI3L1", "MT-CO2", "MT-ND4L", "GADD45B", "MIF", "IFITM3", "CYP1B1"))
dev.copy(png, file = "240709-Mono-40-2000x2000.png", width = 2000, height = 2000)
dev.off() 

FeaturePlot(object = integrated, features = c("THBS1", "PHACTR1", "PYCARD", "MS4A7", "MAP3K8", "DUSP2", "ID1", "MTSS1", "SMPDL3A", "ATF3",   "PTMA", "NPC2", "CCL2", "CTSC",  "EGR1", "IFNGR1", "KLF2", "CPVL",  "MRC1",   "TCF4", "UTRN", "CDKN1A", "PDE4B", "CALR", "ZFP36L1", "PDK4", "IFIT2",  "TPP1",  "BHLHE40", "AZU1",   "BNIP3L", "ID3"))
dev.copy(png, file = "240709-Mono-60-2000x2000.png", width = 2000, height = 2000)
dev.off() 

# MONO: "MS4A4A", "CD69", "S100A4", "NCF1", "EV12B", "C1QA", "MSRB1", "PYCARD", "MAP3K8", "CPVL", "MRC1", "PDK4", "AZU1", "CLEC10A", "PLBD1", "ITGA4", "PTGS2"
FeaturePlot(integrated, features = c("MS4A4A", "CD69", "S100A4", "NCF1", "EVI2B", "C1QA", "MSRB1", "PYCARD", "MAP3K8", "CPVL", "MRC1", "PDK4", "AZU1", "CLEC10A", "PLBD1", "ITGA4", "PTGS2", "THBS1", "PHACTR1", "PYCARD"))
dev.copy(png, file = '240720-MOnocyte-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

# MEGA
FeaturePlot(integrated, features = c("CD2", "IFI6", "HBA2", "F13A1", "HBG2", "IRX3", "TUBB1", "LGALSL", "TUBA4A","PPBP","GP9", "MTURN",  "BEX1","NRGN"))
dev.copy(png, file = '240709-Megakaryocytes-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

# EARLY ERY: "ZFP36", "C1QB", "C1QC","AIF1" ,"ALAS2","HBD", "ATG10",
FeaturePlot(integrated, features = c("ZFP36", "C1QB", "C1QC","AIF1" ,"ALAS2","HBD", "ATG10", "IRX3", "TUBB1", "LGALSL", "TUBA4A","PPBP","GP9", "MTURN",  "BEX1","NRGN",  "CD74","HBA1", "RHCE", "GYPA"))
dev.copy(png, file = '240709-Early-Eryth-2000x2000.png', width = 2000, heigh = 2000)
dev.off()


# MID ERY
FeaturePlot(integrated, features = c( "AC021224.1", "CD74","HBA1", "RHCE", "GYPA", "EEF1A1",  "LDHB", "NCL",   "REXO2", "SYNGR1", "CA1" , "HMGB1",  "RPS12", "STMN1", "H2AFZ", "GAPDH", "CALM1", "RPL22L1", "CDK6", "FAM178B"))
dev.copy(png, file = '240709-Mid-Eryth-20-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(integrated, features = c("SLC25A37", "RPS2", "TMEM141","SLC25A39", "HMBS", "P4HB",  "MGST3", "GLRX", "DSTN", "KCNH2", "LMNA", "NME1", "MT-ATP6", "TXN", "TAGLN2", "NASP", "EPCAM", "PARVB", "CCDC85B", "TSC22D1"))
dev.copy(png, file = '240709-Mid-Eryth-40-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

FeaturePlot(integrated, features = c("MYC", "MT-CO3", "PDIA6", "NOP16",  "MT-ND2",  "CD320", "CNRIP1", "APOC1", "FABP5", "HSPB1", "MT-ND1", "MIF", "CD63",   "JUN",   "GDF15", "CDC20", "MT-ND3", "TMSB10"))
dev.copy(png, file = '240709-Mid-Eryth-60-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

# OUTPUT: 'GYPA', 'KCNH2',  'HBA1','APOC1', 'MYC','CDC20',
FeaturePlot(integrated, features = c( 'GYPA', 'KCNH2',  'HBA1','APOC1', 'MYC','CDC20', "CD2", "IFI6", "HBA2", "F13A1", "HBG2", "IRX3", "TUBB1", "LGALSL", "TUBA4A","PPBP","GP9", "MTURN",  "BEX1","NRGN", 'CDF15'))
dev.copy(png, file = '240709-Mid-Erythroid-2000x2000.png', width = 2000, heigh = 2000)
dev.off()


# Late ERY: "CD14", "MS4A6A", "PLCG2", "LST1", "TYMP", "CSF1R", "HBE1", "CD68"
FeaturePlot(integrated, features = c("CD14", "MS4A6A", "PLCG2", "LST1", "TYMP", "CSF1R", "HBE1", "CD68", 'GYPA', 'KCNH2',  'HBA1','APOC1', 'MYC','CDC20', "CD2", "IFI6", "HBA2", "F13A1", "HBG2"))
dev.copy(png, file = '240709-Late-Erythroid-2000x2000.png', width = 2000, heigh = 2000)
dev.off()


# Erythroid combined
FeaturePlot(integrated, features = c("ZFP36", "C1QB", "C1QC","AIF1" ,"ALAS2","HBD",  'GYPA', 'KCNH2',  'HBA1','APOC1', 'MYC','CDC20',"CD14", "MS4A6A","PLCG2", "LST1", "TYMP", "CSF1R", "HBE1", "CD68"))
dev.copy(png, file = '240709-Erythroid-2000x2000.png', width = 2000, heigh = 2000)
dev.off()


# = = = = = = = = = = = = = = = = = AF 20220 BLOOD MARKERS = = = = = = = = = = = = = = = = = = = = = = = = = = # 240714
# naive (immature, uncommitted progenitors): 'KIT', 'GATA2', + lack of lineage markers
# Megakaryocyte:  'GP9', 'PF4', 'GATA1', 'TAL1', 'FLI1'
# ERYTHROID: 'GYPA', 'KLF1', 'MYC'
# GRANULOCYTE: 'AZU1', 'CEBP-D', 'CEBP-B', 'CEBP-A', 'CEBP-E', 'PRNT3'

FeaturePlot(integrated, features = c('KIT', 'GATA2', 'GP9', 'PF4', 'GATA1', 'TAL1', 'FLI1', 'GYPA', 'KLF1','MYC', 'AZU1', 'CEBPD', 'CEBPB', 'CEBPA', 'CEBPE', 'CEBP'))
dev.copy(png, file = '240712-LINEAGE-2000x2000.png', width = 2000, heigh = 2000)
dev.off()

# ARTERIAL: 'DLL4', 'GJA4', 'GJA5', 'EFNB2', 'SOX17', 'HEY1'
FeaturePlot(integrated, features = c('DLL4', 'GJA4', 'GJA5', 'EFNB2', 'SOX17', 'HEY1'), split.by = 'Condition')
dev.copy(png, file = "24714-arterial-cell-exp-split-by-condition-1095x2000.png", width = 1095, height = 2000)
dev.off()

DotPlot(object = integrated, features = c('DLL4', 'GJA4', 'GJA5', 'EFNB2', 'SOX17', 'HEY1'), group.by = "Condition") +
  ggtitle("ARTERIAL CELLS EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Condition") +
  ylab("Features")

dev.copy(png, file = "24714-arterial-cell-exp-split-by-conditiondotplot-900x700.png", width = 900, height = 700)
dev.off()

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



# 240720
#VENOUS ENDOTHELIU,: 'NR2R2', 'KRT18','FLT4','SMARCA4','EPHB4', 'PCDH12'
#ARTERIAL ENDOTHELIUM: 'DLL4','GJA5','GJA4','EFNB2','SOX17','HEY1'
# ENDOTHELIUM: 'HOPX', 'TEK','
DotPlot(integrated, features = c('HOPX','TEK', # ENDOTHEIUM
                                 'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                                 'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
                                 'CD44','RUNX1','CD82','MYB', # EHT 
                                 'TAL1', 'MEF2C', # MEGA
                                 'AZU1', 'MPO', 'RNASE2', # GRAN
                                 'CD14','PDK4','NCF1', # MONOCYTE
                                 'C1QA','MRC1', # MACROPHAGE
                                 'HLA-DRA','HLA-DRB1'), # IMMUNO-CELLS
        group.by = 'Condition' )+
  ggtitle("GENE EXPRESSION ACROSS CLUSTERS") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  ylab("Conditions")
dev.copy(png, file = "240720-intergrated-expression-distribution-dotplot-CONDITION-500X800.png", width = 500, height = 800)
dev.off()

# BY DAY
DotPlot(integrated, features = c('HOPX', # ENDOTHEIUM
                                 'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                                 'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
                                 'RUNX1','CD82','MYB', # EHT 
                                 'TAL1', 'MEF2C', # MEGA
                                 'AZU1', 'MPO', 'RNASE2', # GRAN
                                 'CD14','PDK4','NCF1', # MONOCYTE
                                 'C1QA','MRC1', # MACROPHAGE
                                 'HLA-DRA','HLA-DRB1'), # IMMUNO-CELLS
        group.by = 'libraryID.1' )+
  ggtitle("GENE EXPRESSION AT DAY 10 AND 13") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  ylab("Conditions")
dev.copy(png, file = "240720-intergrated-expression-distribution-dotplot-DAY-500X800.png", width = 500, height = 800)
dev.off()
