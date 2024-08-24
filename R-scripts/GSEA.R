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
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# Set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")
integrated <- readRDS("240719-integrated-with-cc-data.RDS") 
temp_idents <- Idents(integrated)
DimPlot(integrated)
# 
# prestim_vs_laminar <- FindMarkers(integrated, ident.1 = "Prestimulation", ident.2 = "Laminar")
# prestim_vs_static <- FindMarkers(integrated, ident.1 = "Prestimulation", ident.2 = "Static")
# laminar_vs_static <- FindMarkers(integrated, ident.1 = "Static", ident.2 = "Laminar")
# 
# # Write CSV files
# write.csv(prestim_vs_laminar, "240720-prestim-vs-laminar-findmarkers-for-GSEA.csv", row.names = TRUE)
# write.csv(prestim_vs_static, "240720-prestim-vs-static-findmarkers-for-GSEA.csv", row.names = TRUE)
# write.csv(laminar_vs_static, "240720-laminar-vs-static-findmarkers-for-GSEA.csv", row.names = TRUE)

# Function to create a custom dot plot
custom_dotplot <- function(gsea_results, showCategory = 20, title = "Dotplot", output_file = "dotplot.png") {
  gsea_df <- as.data.frame(gsea_results)
  gsea_df <- gsea_df[1:showCategory,]
  
  # Calculate gene ratio
  gsea_df$GeneRatio <- sapply(gsea_df$GeneRatio, function(x) eval(parse(text = x)))
  
  p <- ggplot(gsea_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = setSize, color = qvalue)) +
    geom_point(alpha = 0.8) +
    scale_color_gradient(low = "red", high = "blue") +
    scale_size(range = c(3, 10)) +
    labs(title = title, x = "Gene Ratio", y = "Gene Set", color = "p.adjust", size = "Gene Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
  
  print(p)
  # Save plot using ggsave
  ggsave(filename = output_file, plot = p, width = 10, height = 8)
  # Save plot using dev.copy
  dev.copy(png, file = output_file, width = 1000, height = 800)
  dev.off()
}

# Define function for GSEA analysis and plotting
run_gsea_and_plot <- function(markers_file, output_prefix) {
  markers <- read.csv(file = markers_file)
  gene_list <- markers$avg_log2FC
  names(gene_list) <- markers$X
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_results <- gseGO(geneList = gene_list, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "SYMBOL", 
                        ont = "BP", 
                        minGSSize = 10, 
                        maxGSSize = 500, 
                        pvalueCutoff = 0.05, eps=0)
  
  custom_dotplot(gsea_results, showCategory = 20, title = paste0(output_prefix, " GSEA Dotplot"), output_file = paste0(output_prefix, "-GSEA-dotplot.png"))
  
  gsea_df_sorted <- as.data.frame(gsea_results)[order(gsea_results$NES, decreasing = TRUE),]
  write.csv(gsea_df_sorted, file = paste0(output_prefix, '-GSEA-sort-by-NES.csv'), row.names = TRUE)
}

# Run GSEA analysis and plot for each comparison
run_gsea_and_plot('ouput(res0.3)/GSEA/240720-prestim-vs-laminar-findmarkers-for-GSEA.csv', "240721-prestim_vs_laminar")
run_gsea_and_plot('ouput(res0.3)/GSEA/240720-prestim-vs-static-findmarkers-for-GSEA.csv', "240721-prestim_vs_static")
run_gsea_and_plot('ouput(res0.3)/GSEA/240720-static-vs-laminar-findmarkers-for-GSEA.csv', "240721-static_vs_laminar")

#======================================================================================================#
static_vs_laminar <- read.csv(file ='ouput(res0.3)/GSEA/240720-static-vs-laminar-findmarkers-for-GSEA.csv')

gene_list_static_vs_laminar<- static_vs_laminar$avg_log2FC
names(gene_list_static_vs_laminar) <- static_vs_laminar$X
gene_list_static_vs_laminar<- sort(gene_list_static_vs_laminar, decreasing = TRUE)
#write.csv(gene_list_static_vs_laminar,file = '240808_gene-list(static_vs_laminar_all_clusters)-GSEA.csv' )
# gene_list_static_vs_laminar<- read.csv(file = '240808_gene-list(static_vs_laminar)-GSEA.csv')

gsea_results <- gseGO(geneList = gene_list_static_vs_laminar, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",   eps = 0)

write.csv(gsea_results, file = 'output-dissertation-figures_240720/240808_static_vs_laminar(all clusters)_GSEA_Results.csv')

# Dotplot of GSEA results
dotplot(gsea_results, showCategory = 10, font.size = 14)
dev.copy(png, file = 'output-dissertation-figures_240720/240808-GSEA-dotplot-static_vs_laminar(all)-show10(up)-660x564.png', width = 660, height = 564)
dev.off()

# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
categorys <- c("isoprenoid biosynthetic process", 
               "nonribosomal peptide biosynthetic process",
               "regulation of phospholipid biosynthetic process", 
               "glutathione biosynthetic process",
               "positive regulation of dendrite extension",
               "glutathione metabolic process",
               "motor neuron apoptotic process",
               "regulation of insulin-like growth factor receptor signaling pathway",
               "cholesterol biosynthetic process",
               "protein localization to cilium")
dotplot(gsea_results, showCategory = categorys, font.size = 14)
dev.copy(png, file = 'output-dissertation-figures_240720/240808-GSEA-dotplot-static_vs_laminar(all)-show10(DOWN)-660x564.png', width = 660, height = 564)
dev.off()


#=======================================================================================
# prestim vs lam
prestim_vs_laminar <- read.csv(file ='ouput(res0.3)/GSEA/240720-prestim-vs-laminar-findmarkers-for-GSEA.csv')
gene_list_prestim_vs_laminar<- prestim_vs_laminar$avg_log2FC
names(gene_list_prestim_vs_laminar) <- prestim_vs_laminar$X
gene_list_prestim_vs_laminar<- sort(gene_list_prestim_vs_laminar, decreasing = TRUE)
#write.csv(gene_list_prestim_vs_laminar,file = '240808_gene-list(prestim_vs_laminar_all_clusters)-GSEA.csv' )
# gene_list_prestim_vs_laminar<- read.csv(file = '240808_gene-list(prestim_vs_laminar)-GSEA.csv')

gsea_results <- gseGO(geneList = gene_list_prestim_vs_laminar, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",   eps = 0)

write.csv(gsea_results, file = 'output-dissertation-figures_240720/240808_prestim_vs_laminar(all clusters)_GSEA_Results.csv')

# Dotplot of GSEA results
dotplot(gsea_results, showCategory = 10, font.size = 14)
dev.copy(png, file = 'output-dissertation-figures_240720/240808-GSEA-dotplot-prestim_vs_laminar(all)-show10(up)-660x564.png', width = 660, height = 564)
dev.off()

# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
categorys <- c("cellular response to cAMP", 
               "glutathione metabolic process",
               "glutathione biosynthetic process", 
               "response to laminar fluid shear stress",
               "TORC2 signaling",
               "blood vessel endothelial cell migration",
               "glycoprotein metabolic process",
               "endothelial cell migration",
               "phospholipid metabolic process",
               "regulation of vasculature development")
dotplot(gsea_results, showCategory = categorys, font.size = 14)
dev.copy(png, file = 'output-dissertation-figures_240720/240808-GSEA-dotplot-prestim_vs_laminar(all)-show10(DOWN)-660x564.png', width = 660, height = 564)
dev.off()

#=======================================================================================
# prestim vs static
prestim_vs_static <- read.csv(file = 'ouput(res0.3)/GSEA/240720-prestim-vs-static-findmarkers-for-GSEA.csv')

gene_list_prestim_vs_static<- prestim_vs_static$avg_log2FC
names(gene_list_prestim_vs_static) <- prestim_vs_static$X
gene_list_prestim_vs_static<- sort(gene_list_prestim_vs_static, decreasing = TRUE)

gsea_results <- gseGO(geneList = gene_list_prestim_vs_static, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "SYMBOL",   eps = 0)

write.csv(gsea_results, file = 'output-dissertation-figures_240720/240808_prestim_vs_static(all clusters)_GSEA_Results.csv')

# Dotplot of GSEA results
pre_cat <- c("ribosomal large subunit biogenesis", 
"negative regulation of sister chromatid segregation",
"nuclear DNA replication", 
"ribosome assembly",
"mitotic spindle assembly checkpoint signaling",
"mitotic spindle checkpoint signaling",
"cell cycle DNA replication",
"mchromosome separation",
"centromere complex assembly",
"mitotic recombination")

dotplot(gsea_results, showCategory = pre_cat, font.size = 14)
dev.copy(png, file = 'output-dissertation-figures_240720/240808-GSEA-dotplot-prestim_vs_static(all)-show10(up)-660x564.png', width = 660, height = 564)
dev.off()

# use `showCategory` to select the displayed terms. It can be a number of a vector of terms.
categorys <- c("cellular defense response", 
               "humoral immune response",
               "antigen processing and presentation", 
               "hydrogen peroxide metabolic process",
               "granulocyte activation",
               "neutrophil activation",
               "regulation of myeloid leukocyte mediated immunity",
               "myeloid cell activation involved in immune response",
               "interleukin-10 production",
               "T cell migration")
dotplot(gsea_results, showCategory = categorys, font.size = 14)
dev.copy(png, file = 'output-dissertation-figures_240720/240808-GSEA-dotplot-prestim_vs_static(all)-show10(DOWN)-660x564.png', width = 660, height = 564)
dev.off()













































# need to convert the gene name into fucking numbers first
gsea_results_KEGG <- gseKEGG(geneList = gene_list_prestim_vs_static, 
                             organism = 'hsa', keyType = 'kegg' 
                             ,eps = 0)
