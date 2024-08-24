# compared laminar and static ECs 
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(glmGamPoi)
library(Seurat)
library(SeuratWrappers)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(monocle3)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(ReactomePA)

# set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")

integrated <- readRDS("RDS/240719-integrated-with-cc-data.RDS") 

Idents(integrated) <- 'seurat_clusters'

# Create a new identity class for Laminar (combining clusters 3 and 5)
integrated$NewIdent <- "Other"  
integrated$NewIdent[WhichCells(integrated, idents = c('3', '5'))] <- "Laminar"

# Create a new identity class for Static (cluster 2)
integrated$NewIdent[WhichCells(integrated, idents = '1')] <- "Static"

# Set the identity of the Seurat object to the new identity class
Idents(integrated) <- "NewIdent"

# Use FindAllMarkers to compare Laminar and Static
Idents(EC_con)
static_vs_laminar <- FindMarkers(EC_con, ident.1 = "Static", ident.2 = "Laminar")

filtered <- static_vs_laminar %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
write.csv(filtered, file = '240808-findmarker_static(1)_EC-vs_laminar(35)_EC.csv', row.names = T)


# = UNO reverse
Idents(EC_con)
# let's use laminar vs static
laminar_vs_static <- FindMarkers(EC_con, ident.1 = "Laminar", ident.2 = "Static")

filtered <- laminar_vs_static %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
# write.csv(filtered, file = 'output-dissertation-figures_240720/EC comparison (static_vs_laminar)/UNO-reverse/240809-findmarker_laminar(35)_EC-vs_static(1)_EC.csv', row.names = T)
laminar_vs_static <- read.csv(file = 'output-dissertation-figures_240720/EC comparison (static_vs_laminar)/UNO-reverse/240809-findmarker_laminar(35)_EC-vs_static(1)_EC.csv')


gene_list_laminar_vs_static <- laminar_vs_static$avg_log2FC
names(gene_list_laminar_vs_static) <- laminar_vs_static$X
gene_list_laminar_vs_static <- sort(gene_list_laminar_vs_static, decreasing = TRUE)
write.csv(gene_list_laminar_vs_static,file = '240809_gene-list(laminar_vs_static)-GSEA.csv' )
# gene_list_laminar_vs_static <- read.csv(file = '240808_gene-list(laminar_vs_static)-GSEA.csv')

#remove duplicate
gene_list_laminar_vs_static <- gene_list_laminar_vs_static[!duplicated(names(gene_list_laminar_vs_static))]

# Set parallel computation to NULL (single-threaded) to avoid connection issues
BiocParallel::register(BiocParallel::SerialParam())

gsea_results <- gseGO(geneList = gene_list_laminar_vs_static, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",   eps = 0)

write.csv(gsea_results, file = '240809_laminar(35)_EC-vs_static(1)_EC_GSEA_Results.csv')

# Dotplot of GSEA results
dotplot(gsea_results, showCategory = 20)
dev.copy(png, file = '240808-GSEA-dotplot-static(1)_EC-vs_laminar(35)_EC-800x1000.png', width = 800, height = 1000)
dev.off()




































#======================================================================================================#
static_vs_laminar <- read.csv(file = '240808-findmarker_static(1)_EC-vs_laminar(35)_EC.csv')

gene_list_static_vs_laminar <- static_vs_laminar$avg_log2FC
names(gene_list_static_vs_laminar) <- static_vs_laminar$X
gene_list_static_vs_laminar <- sort(gene_list_static_vs_laminar, decreasing = TRUE)
write.csv(gene_list_static_vs_laminar,file = '240808_gene-list(static_vs_laminar)-GSEA.csv' )
# gene_list_static_vs_laminar <- read.csv(file = '240808_gene-list(static_vs_laminar)-GSEA.csv')

gsea_results <- gseGO(geneList = gene_list_static_vs_laminar, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",   eps = 0)

write.csv(gsea_results, file = '240808_static(1)_EC-vs_laminar(35)_GSEA_Results.csv')

# Dotplot of GSEA results
dotplot(gsea_results, showCategory = 20)
dev.copy(png, file = '240808-GSEA-dotplot-static_vs_laminar-800x1000.png', width = 800, height = 1000)
dev.off()

# use KEGG by reactome package - fucking hell it kept saying no gene can be mapped
# I even converted the symbol to fucking id yet now it kept saying HTTP 404 bad fucking request
# i will use shinygo instead fuck you all

# Convert gene symbols to Entrez IDs
gene_symbols <- names(gene_list_static_vs_laminar)  # Assuming your geneList is named
gene_df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Ensure the gene list has Entrez IDs
gene_list_entrez <- gene_list_static_vs_laminar[gene_df$SYMBOL]
names(gene_list_entrez) <- gene_df$ENTREZID

gsea_results_KEGG <- gseKEGG(geneList = gene_list_entrez, 
                             organism = 'hsa', 
                             keyType = "ENTREZID", 
                             eps = 0)

# write.csv(gsea_results, file = '240808_static(1)_EC-vs_laminar(35)_GSEA_Results-KEGG.csv')

# Dotplot of GSEA results
dotplot(gsea_results_KEGG, showCategory = 20)
dev.copy(png, file = '240808-GSEA-dotplot-static_vs_laminar-KEGG-800x1000.png', width = 800, height = 1000)
dev.off()









# LAMINAR: "FAM167B", "KLF4", "IGFBP5", 'STMN2','UNC13A','FABP3','KLF2',
# STATIC:

FeaturePlot(integrated, feature = c("AL161629.1", "AC137932.2", "INHBB", "HS3ST2", "DDIT4L", "NCALD", "KLF9", "ZFYVE28", "DEPTOR", "PTHLH", "GGT5", "G0S2", "SNTB1", "SPRY4", "LAMA4", "CAPN8", "LINC02202", "LYVE1", "SLC24A2", "AJAP1", 
                                    "AC012447.1", "DKK1", "ATP8B4", "DGKB", "MS4A4A", "LINC01629", "CGN", "PCLAF", "PCDH18", "BRIP1", "DIO3", "GLP1R", "HIST1H3B", "FAT3", "SMAD9", "AL022328.4", "CMTM5", "DIAPH3", "AL021155.5", "POLQ"))
dev.copy(png, file = "24808-EC-static_laminar-featureplot.png", width = 800, height = 1095)
dev.off()

DimPlot(integrated, split.by = 'Condition', group.by = 'seurat_clusters', label = T)
