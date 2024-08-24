# from blood 2020 script on using GO to find a list of TF + membrane markers
library(Seurat)
library(dplyr)
library(Matrix)
library(tibble)
library(biomaRt)
library(ggplot2)


############################################
# DETERMINATION OF MEMBRANE AND TF MARKERS #
############################################


seuratobject <- readRDS("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/RDS/240721-cluster6_7_pseudotime_included.RDS") 
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/output-dissertation-figures_240720/GO-TF-Membrane-Markers/markers/cluster6+7/")
temp_idents <- Idents(seuratobject)
Idents(seuratobject) <- 'SCT_snn_res.0.2'  

seuratobject <- cluster6_7
Idents(cluster6_7)


# after generating all the csv files =======================================================

# tf: "#8cc2e4", 
# mm: "#af8fd2",
# SELECT ONLY TOP 30
# Define a function to create and save bar plots
plot_and_save_barplot <- function(csv_file, plot_title, output_file, x_limit, label_size, axis_label_size, y_label_size, num_genes = 30) {
  # Read the CSV file
  genelist <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  # Sort the gene list by avg_log2FC in decreasing order
  genelist <- genelist[order(genelist$avg_log2FC, decreasing = TRUE), ]
  
  # Select the top 'num_genes' genes
  genelist <- head(genelist, num_genes)
  
  # Create the bar plot
  png(filename = output_file, width = 700, height = 600)
  barplot(rev(genelist$avg_log2FC), 
          names.arg = rev(genelist$gene), las = 1, horiz = TRUE,
          xlim = x_limit, 
          col =  "#8cc2e4", 
          main = plot_title, 
          cex.names = label_size, 
          cex.axis = y_label_size)
  dev.off()
}

# List of cluster numbers (update this list based on your actual cluster numbers)
cluster_numbers <- c(0,3,4)

# Loop through each cluster number and generate plots
for (clusternumber in cluster_numbers) {
  # Define file paths and plot titles
  tf_csv_file <- paste0("TF_cluster_", clusternumber, ".csv")
  tf_plot_title <- paste0("TOP30 TF Cluster ", clusternumber)
  tf_output_file <- paste0("TOP30_TF_cluster_", clusternumber, ".png")
  
  membrane_csv_file <- paste0("Membrane_cluster_", clusternumber, ".csv")
  membrane_plot_title <- paste0("TOP30 Membrane Marker Cluster ", clusternumber)
  membrane_output_file <- paste0("TOP30_Membrane_cluster_", clusternumber, ".png")
  
  # Plot and save TF markers with x-axis label size of 0.7, y-axis label size of 0.7, and top 30 genes
  plot_and_save_barplot(tf_csv_file, tf_plot_title, tf_output_file, c(0, 5), 0.7, 0.5, 1, 30)
  
  # Plot and save membrane markers with x-axis label size of 0.7, y-axis label size of 0.7, and top 30 genes
  #plot_and_save_barplot(membrane_csv_file, membrane_plot_title, membrane_output_file, c(0, 8), 0.7, 0.5, 1, 30)
}
#  x_limit, label_size, axis_label_size, y_label_size, num_genes
# TF 9 + X LIMT > 8
# MM 679 > 10

# =====================================================

# seuratobject.markers <- FindAllMarkers(seuratobject,
#                                        logfc.threshold = 0.25,
#                                        min.pct = 0.1,
#                                        only.pos = TRUE,
#                                        slot = 'counts')
# 
# for (clusternumber in unique(Idents(seuratobject)))
# {
#   clustermarkers <- seuratobject.markers[seuratobject.markers$cluster == clusternumber & seuratobject.markers$p_val_adj < 0.05,
#                                           c("gene","avg_log2FC", "p_val_adj")]
#   write.csv (clustermarkers, 
#              file = paste0 ("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/output-dissertation-figures_240720/GO-TF-Membrane-Markers/markers/cluster_markers_", clusternumber, ".csv"), 
#              row.names = FALSE)
# }
# rm(clustermarkers)

######## Filtering TF and membrane markers ############
# for (clusternumber in unique(Idents(seuratobject))) 
# {
#   genelist <- read.csv(file = paste0("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/output-dissertation-figures_240720/GO-TF-Membrane-Markers/markers/cluster_markers_", clusternumber, ".csv"),
#                        stringsAsFactors = FALSE)
#   
#   # We are using the Homo Sapiens mart
#   ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#   
#   # GO IDS we are interested in
#   go_id_TF = c('GO:0003700', 'GO:0000130', 'GO:0001071')
#   go_id_membr = c('GO:0005886')
#   #go_id_receptor = 'GO:0098802'
#   filteredgenes <- getBM(attributes = c('hgnc_symbol'),
#                          filters = c('go', 'hgnc_symbol'),
#                          values = list(go_id_TF, genelist$gene),
#                          mart = ensembl, uniqueRows = T)
#   
#   # We can create a subset of values of the original gene list
#   # Find the lines in the original list of genes which contain the result of our query
#   # %in% returns for each position in genelist$gene
#   # TRUE of FALSE if they are or not present in filterdgenes$hgnc-symbol
#   TF.ids <- genelist$gene %in% filteredgenes$hgnc_symbol
#   # Just get those lines
#   genelist.TF <- genelist[TF.ids,]
#   
#   # reassign to genelist.TF the values ordered by logFC
#   genelist.TF <- genelist.TF[order(genelist.TF$avg_log2FC, decreasing = TRUE),]
#   
#   # filter for log FC >=1
#   # genelist.TF <- genelist.TF [(genelist.TF$avg_log2FC >= 1),]
#   
#   write.csv(genelist.TF, file = paste0("TF_cluster_", clusternumber, ".csv"), row.names = FALSE)
#   
#   # barplot ( rev(genelist.TF$avg_log2FC),
#   #           names.arg = rev(genelist.TF$gene),
#   #           las = 1,
#   #           horiz = TRUE,
#   #           xlim = c(0,1.5),
#   #           col = "black",
#   #           main = paste0 ("TF Cluster ", clusternumber),
#   #           cex.names = 1,
#   #           cex.axes = 1)
#   
#   # Save the TF barplot as PNG
#   png(filename = paste0("TF_cluster_", clusternumber, ".png"), width = 800, height = 600)
#   barplot(rev(genelist.TF$avg_log2FC), 
#           names.arg = rev(genelist.TF$gene), 
#           las = 1, 
#           horiz = TRUE,
#           xlim = c(0, 1.5), 
#           col = "#af8fd2", 
#           main = paste0("TF Cluster ", clusternumber), 
#           cex.names = 1, 
#           cex.axes = 1)
#   dev.off()
#   
#   
#   filteredgenes <- getBM(attributes = c('hgnc_symbol'),
#                          filters = c('go', 'hgnc_symbol'),
#                          values = list(go_id_membr, genelist$gene),
#                          mart = ensembl, uniqueRows = T)
#   
#   membr.ids <- genelist$gene %in% filteredgenes$hgnc_symbol
#   # Just get those lines
#   genelist.membr <- genelist[membr.ids,]
#   
#   genelist.membr <- genelist.membr[order(genelist.membr$avg_log2FC, decreasing = TRUE),]
#   
#   write.csv ( genelist.membr, file = paste0("Membrane_cluster_", clusternumber, ".csv"), row.names = FALSE)
#   
#   # barplot ( rev(genelist.membr$avg_log2FC),
#   #           names.arg = rev(genelist.membr$gene),
#   #           las = 1,
#   #           horiz = TRUE,
#   #           xlim = c(0,4),
#   #           col = "#8cc2e4",
#   #           main = paste0("Membrane Marker Cluster ", clusternumber),
#   #           cex.names = 0.3,
#   #           cex.axes = 1)
#   
#   png(filename = paste0("Membrane_cluster_", clusternumber, ".png"), width = 800, height = 600)
#   barplot(rev(genelist.membr$avg_log2FC), 
#           names.arg = rev(genelist.membr$gene), 
#           las = 1, 
#           horiz = TRUE,
#           xlim = c(0, 4), 
#           col = "#8cc2e4", 
#           main = paste0("Membrane Marker Cluster ", clusternumber), 
#           cex.names = 0.3, 
#           cex.axes = 1)
#   dev.off()
# }


# ====================================================
# # generating all markers but sort of got cramped so i didn't save the png in the end
# plot_and_save_barplot <- function(csv_file, plot_title, output_file, x_limit, label_size, axis_label_size, y_label_size) {
#   
#   genelist <- read.csv(csv_file, stringsAsFactors = FALSE)
#   
#   # Create the bar plot
#   png(filename = output_file, width = 1500, height = 4000)
#   barplot(rev(genelist$avg_log2FC), 
#           names.arg = rev(genelist$gene), 
#           las = 1, 
#           horiz = TRUE,
#           xlim = x_limit, 
#           col = "#8cc2e4",
#           main = plot_title, 
#           cex.names = label_size, 
#           cex.axis = y_label_size)
#   dev.off()
# }
# 
# cluster_numbers <- c(0, 1, 2, 3, 4)
# 
# # Loop through each cluster number and generate plots
# for (clusternumber in cluster_numbers) {
#   # Define file paths and plot titles
#   tf_csv_file <- paste0("TF_cluster_", clusternumber, ".csv")
#   tf_plot_title <- paste0("TF Cluster ", clusternumber)
#   tf_output_file <- paste0("TF_cluster_", clusternumber, ".png")
#   
#   membrane_csv_file <- paste0("Membrane_cluster_", clusternumber, ".csv")
#   membrane_plot_title <- paste0("Membrane Marker Cluster ", clusternumber)
#   membrane_output_file <- paste0("Membrane_cluster_", clusternumber, ".png")
#   
#   # Plot and save TF markers with x-axis label size of 0.7 and y-axis label size of 0.7
#   plot_and_save_barplot(tf_csv_file, tf_plot_title, tf_output_file, c(0, 10), 0.8, 0.7, 2)
#   
#   # Plot and save membrane markers with x-axis label size of 0.7 and y-axis label size of 0.7
#   plot_and_save_barplot(membrane_csv_file, membrane_plot_title, membrane_output_file, c(0, 12), 0.8, 0.7, 5)
# }
