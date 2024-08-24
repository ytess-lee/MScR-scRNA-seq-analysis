# ======================== # NOTES # ======================== # 
# cell sorted by CD34+ on day 8 <- origin
# SUBSET STATIC AS TEST SCRIPT
# okay this script will work on static, now i need figure out how to use Day8 as the starting point :)

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

# set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")

# using new RDS file from Anto - 240614
integrated <- readRDS("RDS/240615_pleasework_Mechanobiology_Harmonised_2reps.RDS") 

# ======================== # ASSIGN CLUSTERS BY CONDITIONS # ======================== # 
# temporarly store the cluster idents to be able to re-use them later
temp_idents <- integrated@active.ident

# assign the culture condition to the active identitites for the DEG analysis
Idents(integrated) <- "Condition" 

# Reorder the levels of the condition factor
integrated$Condition <- factor(integrated$Condition, levels = c("Prestimulation", "Static", "Laminar"))
levels(integrated)
View(integrated@meta.data)
integrated <- FindNeighbors(integrated, dims = 1:13)
integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.5)
DimPlot(integrated, reduction = 'umap')


# = = = = = = = = = = = = = = = = =  = = = = = = = = = = = = =  = = = = = = = = = = Annotation = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # 
# Define new cluster IDs
# Detailed annotation
DimPlot(integrated, group.by = 'SCT_snn_res.0.3', label = T)
new.cluster.ids.endo <- c('Endo 1 - enriched in DNA replication',
                          'Endo 2 - enriched in glycan degradation',
                          'Endo 3 - enriched in oxidative phosphorylation',
                          'Endo 4 - enriched in metabolic pathways',
                          'Endo 5 - enriched in ECM-receptor interaction',
                          'Endo 6 - enriched in Notch signalling pathway',
                          'Haematopoietic lineage',
                          'EHT',
                          'Megakaryocyte-committed',
                          'Erythroid-committed')

# Ensure the new cluster IDs have names that match current cluster levels
current_clusters_endo <- levels(integrated@meta.data$SCT_snn_res.0.3)
if (length(current_clusters_endo) != length(new.cluster.ids.endo)) {
  stop("Number of new cluster IDs does not match the number of current clusters")
}
names(new.cluster.ids.endo) <- current_clusters_endo

# Map the new cluster IDs to the existing SCT_snn_res.0.5 column
integrated@meta.data$SCT_snn_res.0.3 <- new.cluster.ids.endo[as.character(integrated@meta.data$SCT_snn_res.0.3)]

# Convert to factor with proper levels
integrated@meta.data$SCT_snn_res.0.3 <- factor(integrated@meta.data$SCT_snn_res.0.3, levels = unique(new.cluster.ids.endo))

DimPlot(integrated, reduction = "umap",group.by = 'SCT_snn_res.0.3', pt.size = 0.5, label.size=5) + NoAxes() + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240722-umap-annotated-clusters-600x502.png", width = 600, height = 502)
dev.off()

# Simplified annotation
DimPlot(integrated, group.by = 'SCT_snn_res.0.3', label = T)
new.cluster.ids.endo <- c('Endo 1',
                          'Endo 2',
                          'Endo 3',
                          'Endo 4',
                          'Endo 5',
                          'Endo 6',
                          'Haematopoietic lineage',
                          'EHT',
                          'Megakaryocyte-committed',
                          'Erythroid-committed')


# Ensure the new cluster IDs have names that match current cluster levels
current_clusters_endo <- levels(integrated@meta.data$SCT_snn_res.0.3)
if (length(current_clusters_endo) != length(new.cluster.ids.endo)) {
  stop("Number of new cluster IDs does not match the number of current clusters")
}
names(new.cluster.ids.endo) <- current_clusters_endo

# Map the new cluster IDs to the existing SCT_snn_res.0.5 column
integrated@meta.data$SCT_snn_res.0.3 <- new.cluster.ids.endo[as.character(integrated@meta.data$SCT_snn_res.0.3)]

# Convert to factor with proper levels
integrated@meta.data$SCT_snn_res.0.3 <- factor(integrated@meta.data$SCT_snn_res.0.3, levels = unique(new.cluster.ids.endo))

DimPlot(integrated, reduction = "umap",group.by = 'SCT_snn_res.0.3', pt.size = 0.5, label.size=5) + NoAxes() + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240722-umap-annotated-clusters(Simplified)-478x433.png", width = 478, height = 433)
dev.off()

DimPlot(integrated, reduction = "umap",group.by = 'SCT_snn_res.0.3', pt.size = 0.5, label.size=3.5, label = T) + NoAxes() + ggtitle(NULL) + NoLegend()
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240722-umap-annotated-clusters(Simplified)-label-TRUE-478x502.png", width = 478, height = 502)
dev.off()


# SUBSET CLUSTER 6 + 7 Annotation #
cluster6_7 <- readRDS(file = 'RDS/240714-cluster6+7(EHT+Blood)-subset.RDS')
cluster6_7 <- FindClusters(cluster6_7, dims= 1:15, resolution = 0.2)
DimPlot(cluster6_7, group.by = 'SCT_snn_res.0.2', label = T)
new.cluster.ids <- c('Endothelial cells primed for EHT', 
                     'Myeloid-committed',
                     'Immuno-committed cells',
                     'EHT', 
                     'Endothelial cells')

# Ensure the new cluster IDs have names that match current cluster levels
current_clusters <- levels(cluster6_7@meta.data$seurat_clusters)
if (length(current_clusters) != length(new.cluster.ids)) {
  stop("Number of new cluster IDs does not match the number of current clusters")
}
names(new.cluster.ids) <- current_clusters

# Map the new cluster IDs to the existing SCT_snn_res.0.5 column
cluster6_7@meta.data$seurat_clusters <- new.cluster.ids[as.character(cluster6_7@meta.data$seurat_clusters)]

# Convert to factor with proper levels
#cluster6_7@meta.data$SCT_snn_res.0.2 <- factor(cluster6_7@meta.data$SCT_snn_res.0.5, levels = unique(new.cluster.ids))
DimPlot(cluster6_7, reduction = 'umap', group.by = "seurat_clusters", pt.size = 0.5, label.size=3.5) + NoAxes() + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240722-umap-annotated-EHT-to-lineage-478X502.png", width = 478, height = 502)
dev.off()

# ECs broad annotation
DimPlot(integrated, reduction = "umap",group.by = 'seurat_clusters', pt.size = 0.5, label.size=3.5) + NoAxes() + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240720-umap-broad-ECs-annotated-478X502.png", width = 478, height = 502)
dev.off()


# = = = = = = = = = = = = = = = = = = = = = = = = Trajectory analysis =  = = = = = = = = = = = = = # 
####### PSEUDOTEMPORAL ORDERING ######
DimPlot(integrated, group.by = 'SCT_snn_res.0.3')
integrated$Condition <- factor(integrated$Condition, levels = c("Prestimulation", "Static", "Laminar"))
Idents(integrated) <- "Condition"

prestim_vs_static <- subset(integrated, idents = c('Prestimulation','Static'))
saveRDS(prestim_vs_static, file = "240721-prestim_vs_static-subset-for-trajectory.RDS") # you can't assign the directory of where it's saved to using this func

# Prestim vs Static
monocle_object <- as.cell_data_set(prestim_vs_static) # Warning: Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object
monocle_object <- cluster_cells(cds = monocle_object, reduction_method = "UMAP")
monocle_object <- learn_graph(monocle_object, use_partition = FALSE)
DimPlot(prestim_vs_static, split.by = "libraryID.1") + NoAxes() 
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_static_split-day-478X502.png", width = 478, height = 502)
dev.off()

# here I choose the further up right point as the starting node
monocle_object <- order_cells(monocle_object, reduction_method = "UMAP")

#monocle_object <- reduce_dimension(monocle_object)
plot_cells(monocle_object, 
           label_groups_by_cluster=TRUE,  
           color_cells_by = "SCT_snn_res.0.3",
           group_label_size = 4)+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_static_plot_cells(AF + new starting pt)-478X502.png", width = 478, height = 502)
dev.off()

plot_cells(monocle_object,
           color_cells_by = "SCT_snn_res.0.3",
           label_groups_by_cluster = T, 
           group_label_size = 4) +
  scale_color_manual(values = c("#8cc2e4", "#90c47d", "#af8fd2", "#e48484", "#e4ad84", "#e4e484", "#8ce1e1", "#e48398", "#c2c2c2", "#deb698"))+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_static_plot_cells(new start pt)-478X502.png", width = 478, height = 502)
dev.off()


plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           show_trajectory_graph = FALSE,
           cell_size = 0,
           cell_stroke = 2,
           alpha = 0.5) + NoAxes()
dev.copy(png, file = "dissertation-figures_240720-prestim_vs_static_plot_cells_pseudotime(AF +(new start pt))-478X502.png", width = 478, height = 502)
dev.off()


plot_cells(monocle_object,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 7)+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_static_plot_cells_pseudotime(new start pt)-478X502.png", width = 478, height = 502)
dev.off()

prestim_vs_static$pseudotime <- pseudotime(monocle_object)
monocle_object@clusters$seurat <- prestim_vs_static@meta.data$SCT_snn_res.0.3

# function to plot features voer pseudotime
# all coloured
plot_features_over_pseudotime <- function(seurat_obj, features,
                                          pseudotime = NULL, loess_span = 0.5, show_se = FALSE) {
  
  # Plots the expression of a set of features over pseudotime.
  # :param seurat_obj: Seurat object.
  # :param features: Vector of feature names.
  # :param pseudotime: Vector of pseudotime values. If NULL, the pseudotime values are taken from the pseudotime metadata of the Seurat object.
  # :param loess_span: Span parameter for the LOESS smoothing function.
  # :param show_se: Whether to show the error bands around the LOESS curve.
  # :return: ggplot object.
  
  if (is.null(pseudotime)) {
    # Check if pseudotime metadata exists
    if (!"pseudotime" %in% names(seurat_obj@meta.data)) {
      stop("No pseudotime metadata found in Seurat object. Please provide pseudotime values using the pseudotime argument.")
    }
    
    pseudotime <- seurat_obj$pseudotime
  }
  
  df <- data.frame(Pseudotime = pseudotime,
                   Genes = t(as.matrix(GetAssayData(seurat_obj)[features,])))
  print(head(df))
  
  df %>% 
    # Transform to long format
    pivot_longer(-Pseudotime, names_to = "Gene", values_to = "Expression") %>%
    # Remove "genes." from the gene names
    mutate(Gene = gsub("Genes.", "", Gene)) %>%
    ggplot(aes(x = Pseudotime, y = Expression, color = Gene)) +
    geom_smooth(method = "loess", se = show_se, span = loess_span) +
    theme_bw() +
    theme(text = element_text(size = 18)) +
    labs(x = "Pseudotime", y = "Expression")
}


plot_features_over_pseudotime(prestim_vs_static, features = c('CDK4','CDK6','CCNE2',
                                                              'CDK2','CCND3', 'CCNA2',
                                                              'CDK1', 'CCNB1'), show_se = T)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_static_pseudotime-func-test-cc_markers-800X400.png", width = 800, height = 400)
dev.off()


# plot_features_over_pseudotime
# simplified colour coding for better HE/ HPc visualisation plot
plot_features_over_pseudotime <- function(seurat_obj, features,
                                          pseudotime = NULL, loess_span = 0.5, show_se = FALSE) {
  # Plots the expression of a set of features over pseudotime.
  # :param seurat_obj: Seurat object.
  # :param features: Vector of feature names.
  # :param pseudotime: Vector of pseudotime values. If NULL, the pseudotime values are taken from the pseudotime metadata of the Seurat object.
  # :param loess_span: Span parameter for the LOESS smoothing function.
  # :param show_se: Whether to show the error bands around the LOESS curve.
  # :return: ggplot object.
  
  if (is.null(pseudotime)) {
    # Check if pseudotime metadata exists
    if (!"pseudotime" %in% names(seurat_obj@meta.data)) {
      stop("No pseudotime metadata found in Seurat object. Please provide pseudotime values using the pseudotime argument.")
    }
    
    pseudotime <- seurat_obj@meta.data$pseudotime
  }
  
  df <- data.frame(Pseudotime = pseudotime,
                   Genes = t(as.matrix(GetAssayData(seurat_obj)[features,])))
  
  plot_data <- df %>% 
    pivot_longer(-Pseudotime, names_to = "Gene", values_to = "Expression") %>%
    mutate(Gene = gsub("Genes.", "", Gene))
  
  # Assign colours
  plot_data <- plot_data %>%
    mutate(Color = ifelse(Gene %in% c("DLL4", "CDH5"), "#e48398", "#8cc2e4"))
  
  ggplot(plot_data, aes(x = Pseudotime, y = Expression, color = Gene, group = Gene)) +
    geom_smooth(method = "loess", se = show_se, span = loess_span) +
    scale_color_manual(values = c( "DLL4" = "#e48398", "RUNX1" = "#8cc2e4",
                                   "CDH5" = "#e48398", "CD44" = "#8cc2e4", "GATA2" = "#8cc2e4",
                                   "SPINK2" = "#8cc2e4", "SMARCA2" = "#8cc2e4")) +
    theme_bw() +
    theme(text = element_text(size = 18), 
          axis.text.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)) +
    labs(x = "Pseudotime", y = "Expression", color = "Gene")
  # 
  # theme(axis.text.x = element_text(size = 12),  # Change font size for x-axis
  #       axis.text.y = element_text(size = 12),  # Change font size for y-axis
  #       axis.title.x = element_text(size = 14), # Change font size for x-axis title
  #       axis.title.y = element_text(size = 14)) 
  # 
  # 
}

he <- c('CDH5', 'DLL4') 
hpc <- c('CD44', 'GATA2', 'RUNX1', 'SPINK2')


# Call the function with your specific parameters
plot1 <- plot_features_over_pseudotime(cluster6_7, features = c(hpc), show_se = TRUE)
ggsave(filename = "output-dissertation-figures_240720/Trajectory/EHT+Blood/240813-cluster6_7_pseudotime-hpc-only-font.png", plot = plot1, width = 10, height = 3, dpi = 100)

plot2 <- plot_features_over_pseudotime(cluster6_7, features = c(he, hpc), show_se = TRUE)
ggsave(filename = "output-dissertation-figures_240720/Trajectory/EHT+Blood/240812-cluster6_7_pseudotime-func-test-cc_markers-800X300.png", plot = plot2, width = 10, height = 3, dpi = 100)
