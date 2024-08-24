# this script is to generate figures used in ppt and dissertation
# + NoAxes() + ggtitle(NULL)

library(dplyr)
library(tidyverse)
library(patchwork)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggridges)
library(ggrepel)
library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)

Sys.setenv(LANG = "en")
# library(SeuratWrappers)
# library(monocle3)
# library(glmGamPoi)
# set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")


integrated <- readRDS("RDS/240719-integrated-with-cc-data.RDS") 

Idents(integrated) <- 'seurat_clusters'
EC <- subset(integrated, ident = c('1','3','5'))

DimPlot(EC, group.by = 'seurat_clusters')
Idents(EC) <- 'Condition'
Idents(EC)

EC_con <- subset(EC, ident = c('Static','Laminar'))
DimPlot(cluster6_7, group.by = 'seurat_clusters') + NoAxes()
notch_target <- c('DLL4','NOTCH1','NOTCH4','HES1','HEY1','HEY2','EFNB2','CCND1')
#FeaturePlot(EC_con, features = c('JAG1', 'DLL4'))

VlnPlot(EC_con, features = c('KLF2','KLF4','NOS3','SMAD6','SMAD7'), 
        group.by = 'seurat_clusters',
        add.noise = F, pt.size = 0)

integrated <- FindNeighbors(integrated, dims = 1:13)

# = = = = = = = = = = = = = = = = =  = = = = = = = = = = = = = = = = = = Resolution =  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # 

ElbowPlot(integrated, ndims = 50)+
  theme(axis.text.x = element_text(size = 16),  # Change font size for x-axis
        axis.text.y = element_text(size = 16),  # Change font size for y-axis
        axis.title.x = element_text(size = 16), # Change font size for x-axis title
        axis.title.y = element_text(size = 16)) 

integrated <- FindClusters(integrated, dims= 1:13, resolution = 0.3)
temp_idents <- Idents(integrated)

DimPlot(integrated, group.by = 'seurat_clusters', label = T) +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240720-umap-seurat-clusters(13)-478x433.png", width = 478, height = 433)
dev.off()

DimPlot(integrated, group.by = 'Condition') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()+ ggtitle(NULL)
dev.copy(png, file = "240720-umap-seurat-clusters(by condition)-478x433.png", width = 478, height = 433)
dev.off()


# ========================================================== one-time-stay
Idents(integrated)
DimPlot(integrated, split.by = 'Condition') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/240721-umap-seurat-clusters(split-by-Condition)-integrated-478x300.png", width = 478, height = 300)
dev.off()

DimPlot(integrated, split.by = 'libraryID.1') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/240721-umap-seurat-clusters(split-by-Day)-integrated-478x300.png", width = 478, height = 300)
dev.off()

DimPlot(integrated, split.by = 'dataset') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()  + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/240721-umap-seurat-clusters(split-by-dataset)-integrated-478x300.png", width = 478, height = 300)
dev.off()

DimPlot(integrated, group.by = 'Condition') +
  xlab(expression(UMAP[1])) +
  ylab(expression(UMAP[2])) + NoAxes()+ ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/240721-umap-seurat-clusters(by condition)-integrated-478x300.png", width = 478, height = 300)
dev.off()




# = = = = = = = = = = = = = = = = = = Feature visualisation = = = = = = = = = = = = = = = = = = = = = = = # 
#VASCULAR ENDOTHELIU,: 'NR2R2', 'KRT18','FLT4','SMARCA4','EPHB4', 'PCDH12'
#ARTERIAL ENDOTHELIUM: 'DLL4','GJA5','GJA4','EFNB2','SOX17','HEY1'
# ENDOTHELIUM: 'HOPX', 'TEK','
DotPlot(integrated, features = c('HOPX','TEK', # ENDOTHEIUM
                                 'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                                 'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
                                 'CD44','RUNX1','CD82','MYB', # EHT 
                                 'TAL1', 'PF4', # MEGA
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

FeaturePlot(cluster6_7, features = c('HOPX','TEK', # ENDOTHEIUM
                                     'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                                     'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
                                     'CD44','RUNX1','CD82','MYB', # EHT 
                                     'TAL1', 'MEF2C', # MEGA
                                     'AZU1', 'MPO', 'RNASE2', # GRAN
                                     'CD14','PDK4','NCF1', # MONOCYTE
                                     'C1QA','MRC1', # MACROPHAGE
                                     'HLA-DRA','HLA-DRB1')) +NoAxes()# IMMUNO-CELLS
dev.copy(png, file = "24721-Endo-EHT-BC-featureplot.png", width = 800, height = 1095)
dev.off()


# BY DAY
DotPlot(cluster6_7, features = c('HOPX', # ENDOTHEIUM
                                 'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                                 'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
                                 'RUNX1','CD82','MYB', # EHT 
                                 'TAL1', 'PF4', # MEGA
                                 'AZU1', 'MPO', 'RNASE2', # GRAN
                                 'CD14','PDK4','NCF1', # MONOCYTE
                                 'C1QA','MRC1', # MACROPHAGE
                                 'HLA-DRA','HLA-DRB1'), # IMMUNO-CELLS
        group.by = 'Condition' )+
  ggtitle("GENE EXPRESSION AT DAY 10 AND 13") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  ylab("Conditions")
dev.copy(png, file = "240720-intergrated-expression-distribution-dotplot-DAY-500X800.png", width = 500, height = 800)
dev.off()

# = = = = = = = = = = = = = = = = =  = = = = = = = = = = = = =  = = = = = = = = = Cell typing =  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # 

# ======================== # SUBSET CLUSTER 6 + 7 # ======================== # 
cluster6_7 <- readRDS(file = 'RDS/240714-cluster6+7(EHT+Blood)-subset.RDS')
cluster6_7 <- FindClusters(cluster6_7, dims= 1:15, resolution = 0.2)

DimPlot(cluster6_7, reduction = "umap") + NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/cell-typing/240721-cluster6+7-subset-478x300.png", width = 478, height = 300)
dev.off()

FeaturePlot(cluster6_7, features = c('TEK','PCDH12', # venous
                                     'SOX17','HEY1','DLL4', # arterial endothelium
                                     'CD44','RUNX1','CD82','MYB', # EHT
                                     'TAL1', # mega
                                     'AZU1', 'MPO', 'RNASE2', # gran
                                     'CD14','PDK4','NCF1', # MONOCYTE
                                     'C1QA','MRC1', # macrophage
                                     'HLA-DRA','HLA-DRB1')) # immuno-committed
dev.copy(png, file = 'output-dissertation-figures_240720/cell-typing/2240721-Cluster6+7-cell-typing-800x800.png', width = 800, heigh = 800)
dev.off()

FeaturePlot(cluster6_7, features = c('CD34','KDR','CDH5'))

# ====
DotPlot(cluster6_7, features = c('TEK','PCDH12', # venous
                                 'SOX17','HEY1','DLL4', # arterial endothelium
                                 'CD44','RUNX1','CD82','MYB', # EHT
                                 'TAL1', # mega
                                 'AZU1', 'MPO', 'RNASE2', # gran
                                 'CD14','PDK4','NCF1', # MONOCYTE
                                 'C1QA','MRC1', # macrophage
                                 'HLA-DRA','HLA-DRB1'), group.by = "seurat_clusters") +
  ggtitle("CELL MARKER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Features") +
  ylab("Sub-clusters")
dev.copy(png, file = "output-dissertation-figures_240720/cell-typing/2240721-CLUSTER6+7-subset-expression-distribution-dotplot-500x650.png", width = 500, height = 650)
dev.off()

# ======================== # SUBSET CLUSTER 8+9 # ======================== # 
Idents(integrated) <- temp_idents
cluster8_9 <- subset(integrated, idents= c("8", "9"))
cluster8_9 <- readRDS(file = 'RDS/240709-cluster-8+9-subset.RDS') # fuck these are flipped
saveRDS(cluster8_9, file = 'RDS/240709-cluster-8+9-subset.RDS')
DimPlot(cluster8_9) + NoAxes()

FeaturePlot(cluster8_9, features = c( 'GP9', 'PF4', 'GATA1', 'TAL1', 'FLI1', 'GYPA', 'KLF1','MYC', 'HBZ', 'HBM'))
dev.copy(png, file = "output-dissertation-figures_240720/cell-typing/2240721-CLUSTER8+9-subset-ery+mega-650x500.png", width = 650, height = 500)
dev.off()

DotPlot(cluster8_9, features = c('GP9', 'PF4', 'GATA1', 'TAL1', 'FLI1',
                                        'GYPA', 'KLF1','MYC', 'HBZ', 'HBM'), 
        group.by = "seurat_clusters") +
  ggtitle("CLUSTER 8, 9 MARKER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Features") +
  ylab("Clusters")
dev.copy(png, file = "output-dissertation-figures_240720/cell-typing/2240721-CLUSTER8+9-subset-ery+mega-dotplot-420X500.png", width = 420, height = 500)
dev.off()

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

# simplified annotation + SPLITY BY conditions
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

DimPlot(integrated, reduction = "umap",group.by = 'SCT_snn_res.0.3', pt.size = 0.5, label.size=3.5, split.by = 'libraryID.1') + NoAxes() + ggtitle(NULL)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240720-umap-annotated-clusters(Simplified+SPLITY-BY-Day)-478X502.png", width = 478, height = 502)
dev.off()

DimPlot(integrated, reduction = "umap",split.by = 'SCT_snn_res.0.2', ncol = 9) + NoAxes() + ggtitle(NULL)

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

# = = = = = = = = = = = = = = = = =  = = = = = = = = = = = = =  = = = = = = = = = = = = = = Trajectory analysis =  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # 

# AF's script

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
saveRDS(prestim_vs_static, file = '240721-prestim_vs_static_pseudotime_included(new start pt).RDS')

# prstim_vs_static
prestim_vs_static <- readRDS(file = '240721-prestim_vs_static_pseudotime_included.RDS')

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

# this is just to check out if the func work with more features, i think i'll need to look for something more specific for this
# i guess i can subset the EHT cells out 
plot_features_over_pseudotime(prestim_vs_static, features = c('CDK4','CDK6','CCNE2',
                                                              'CDK2','CCND3', 'CCNA2',
                                                              'CDK1', 'CCNB1'), show_se = T)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_static_pseudotime-func-test-cc_markers-800X400.png", width = 800, height = 400)
dev.off()



# god knows what it does
pseudotime_df <- data.frame(cell_cycle = prestim_vs_static$CellCycle,
                            pseudotime = prestim_vs_static$pseudotime,
                            RUNX1 = GetAssayData(prestim_vs_static)["RUNX1",], 
                            CDH5 = GetAssayData(prestim_vs_static)["CDH5",],
                            CD44 = GetAssayData(prestim_vs_static)["CD44",],
                            GATA3 = GetAssayData(prestim_vs_static)["GATA3",],
                            GATA2 = GetAssayData(prestim_vs_static)["GATA2",],
                            SPINK2 = GetAssayData(prestim_vs_static)["SPINK2",])

ggplot(pseudotime_df, aes(x = pseudotime, y = CD44)) +
  geom_jitter(aes(col = cell_cycle), size = 0.2, height = 0.01) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  xlab("Pseudotime") +
  ylab("CD44") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

# session restarted, will need to run the monocle_object again for this
plot_genes_in_pseudotime(monocle_object["RUNX1", ],
                         color_cells_by= "CellCycle",
                         min_expr=0.5)




# Prestim vs laminar
prestim_vs_laminar <- subset(integrated, idents = c('Prestimulation','Laminar'))
saveRDS(prestim_vs_laminar, file = "240721-prestim_vs_laminar-subset-for-trajectory.RDS") # haven't saved it yet dw

# Prestim vs Static
monocle_object <- as.cell_data_set(prestim_vs_laminar) # Warning: Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object
monocle_object <- cluster_cells(cds = monocle_object, reduction_method = "UMAP")
monocle_object <- learn_graph(monocle_object, use_partition = FALSE)
DimPlot(prestim_vs_laminar, split.by = "libraryID.1") + NoAxes() 
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_laminar_split-day-478X502.png", width = 478, height = 502)
dev.off()

# here I choose the further up right point as the starting node
monocle_object <- order_cells(monocle_object, reduction_method = "UMAP")

#monocle_object <- reduce_dimension(monocle_object)
plot_cells(monocle_object, 
           label_groups_by_cluster=TRUE,  
           color_cells_by = "SCT_snn_res.0.3",
           group_label_size = 4)+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_laminar_plot_cells(AF)-478X502.png", width = 478, height = 502)
dev.off()

plot_cells(monocle_object,
           color_cells_by = "SCT_snn_res.0.3",
           label_groups_by_cluster = T, 
           group_label_size = 4) +
  scale_color_manual(values = c("#8cc2e4", "#90c47d", "#af8fd2", "#e48484", "#e4ad84", "#e4e484", "#8ce1e1", "#e48398", "#c2c2c2", "#deb698"))+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_laminar_plot_cells-478X502.png", width = 478, height = 502)
dev.off()


plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           show_trajectory_graph = FALSE,
           cell_size = 0,
           cell_stroke = 2,
           alpha = 0.5) + NoAxes()
dev.copy(png, file = "dissertation-figures_240720-prestim_vs_laminar_plot_cells_pseudotime(AF)-478X502.png", width = 478, height = 502)
dev.off()


plot_cells(monocle_object,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 7)+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240721-prestim_vs_laminar_plot_cells_pseudotime-478X502.png", width = 478, height = 502)
dev.off()

prestim_vs_laminar$pseudotime <- pseudotime(monocle_object)
monocle_object@clusters$seurat <- prestim_vs_laminar@meta.data$SCT_snn_res.0.3
saveRDS(prestim_vs_laminar, file = '240721-prestim_vs_laminar_pseudotime_included.RDS')




# traj for EHT + blood ===============================================#
monocle_object <- as.cell_data_set(cluster6_7) # Warning: Monocle 3 trajectories require cluster partitions, which Seurat does not calculate. Please run 'cluster_cells' on your cell_data_set object
monocle_object <- cluster_cells(cds = monocle_object, reduction_method = "UMAP")
monocle_object <- learn_graph(monocle_object, use_partition = FALSE)
DimPlot(cluster6_7, split.by = "libraryID.1") + NoAxes() 
dev.copy(png, file = "output-dissertation-figures_240720/Trajectory/EHT+Blood/240722-cluster6_7_split-day-478X502.png", width = 478, height = 502)
dev.off()

# here I choose the further up right point as the starting node
monocle_object <- order_cells(monocle_object, reduction_method = "UMAP")

#monocle_object <- reduce_dimension(monocle_object)
plot_cells(monocle_object, 
           label_groups_by_cluster=F,  
           color_cells_by = "SCT_snn_res.0.2",
           group_label_size = 4)+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Trajectory/EHT+Blood/240722-cluster6_7_plot_cells(AF)-478X502.png", width = 478, height = 502)
dev.off()

plot_cells(monocle_object,
           color_cells_by = "SCT_snn_res.0.2",
           label_groups_by_cluster = F, 
           group_label_size = 4) +
  scale_color_manual(values = c("#8cc2e4", "#90c47d", "#af8fd2", "#e48484", "#e4ad84"))+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Trajectory/EHT+Blood/240722-cluster6_7_plot_cells-478X502.png", width = 478, height = 502)
dev.off()


plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           show_trajectory_graph = FALSE,
           cell_size = 0,
           cell_stroke = 2,
           alpha = 0.5) + NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Trajectory/EHT+Blood/240722-cluster6_7_plot_cells_pseudotime(AF)-478X502.png", width = 478, height = 502)
dev.off()


plot_cells(monocle_object,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 7)+ NoAxes()
dev.copy(png, file = "output-dissertation-figures_240720/Trajectory/EHT+Blood/240722-cluster6_7_plot_cells_pseudotime-478X502.png", width = 478, height = 502)
dev.off()

cluster6_7$pseudotime <- pseudotime(monocle_object)
monocle_object@clusters$seurat <- cluster6_7@meta.data$SCT_snn_res.0.2
saveRDS(cluster6_7, file = 'RDS/240714-cluster6+7(EHT+Blood)-subset.RDS')
seurat_obj <- readRDS(file = 'RDS/240714-cluster6+7(EHT+Blood)-subset.RDS')
# Featuress over pseudotime oringial script from AF
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
    theme(text = element_text(size = 22),
          axis.text.x = element_text(size = 22),
          axis.title.y = element_text(size = 22)) +
    labs(x = "Pseudotime", y = "Expression")
}
# 
# 
plot_features_over_pseudotime(cluster6_7, features = c(he, hpc), show_se = T)
dev.copy(png, file = "output-dissertation-figures_240720/Annotated/240814-cluster6_7_pseudotime-func-test-cc_markers-1000x300.png", width = 1000, height = 300)
dev.off()
# 
# 
# 
# 
# plot_features_over_pseudotime(cluster6_7, features = c('RUNX1','CDH5','DLL4','MEIS1','CD44','GATA3','GATA2','SPINK2'), show_se = T)
# dev.copy(png, file = "output-dissertation-figures_240720/Trajectory/EHT+Blood/240721-cluster6_7_pseudotime-func-test-cc_markers(AF script)-800X300.png", width = 1000, height = 300)
# dev.off()

# session restarted, will need to run the monocle_object again for this
plot_genes_in_pseudotime(monocle_object["RUNX1", ],
                         label_by_short_name = FALSE,
                         min_expr=0.5)

cluster6_7 <- seurat_obj
# oh this will be fun ========================================
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

# all colours
p3lot_features_over_pseudotime <- function(seurat_obj, features,
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
    theme(text = element_text(size = 22),
          axis.text.x = element_text(size = 22),
          axis.title.y = element_text(size = 22)) +
    labs(x = "Pseudotime", y = "Expression")
}

p3lot_features_over_pseudotime(cluster6_7, features = c(he, hpc), show_se = T)
dev.copy(png, file = 'output-dissertation-figures_240720/240812-integrated-cc-analysis-group-condition-split-cc-font.png', width = 1200, height = 400)
dev.off()


# = = = = = = = = = = = = = = = = =  = = = = = = = = = = = = = = = = = = Cell cycle = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # 
DimPlot(integrated, split.by = 'Condition', group.by = 'CellCycle', cols = c('#F8766D','#00EE17','#619CFF')) + NoAxes + ggtitle(NULL)
dev.copy(png, file = 'output-dissertation-figures_240720/Trajectory/EHT+Blood/240812-he-hpc-pseudotime-all.png', width = 1200, height = 400)
dev.off()

# = = = = = = = = = = = = = = = = =  = = = = = = = = = = = = =  = = = = = = = = = = = = = = X analysis =  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # 
FeaturePlot(cluster6_7, features ='CD84') + NoAxes()
a

FeaturePlot(integrated, features = c('HOPX','TEK', # ENDOTHEIUM
                                     'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                                     'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
                                     'CD44','RUNX1','CD82','MYB', # EHT 
                                     'TAL1', 'MEF2C', # MEGA
                                     'AZU1', 'MPO', 'RNASE2', # GRAN
                                     'CD14','PDK4','NCF1', # MONOCYTE
                                     'C1QA','MRC1', # MACROPHAGE
                                     'HLA-DRA','HLA-DRB1')) +NoAxes()
DimPlot(cluster6_7)

# Arterial = 'EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4'
# venous = 'NR2F2','EPHB4', 'CDH5', 'ESM1', 'FLT4', 'SOX18', 'SMARCA4'

FeaturePlot(integrated, features = c('EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4', # arterial
                                     'NR2F2','EPHB4', 'CDH5', 'ESM1', 'FLT4', 'SOX18', 'SMARCA4')) # venous))
dev.copy(png, file = "output-dissertation-figures_240720/240722-integrated-arterial-venous-endo-markers-800x1095.png", width = 800, height = 1095)
dev.off()

# Idents(integrated) <- 'SCT_snn_res.0.3'
# cluster3 <- subset(integrated, idents = '3')
# cluster5 <- subset(integrated, idents = '5')

# FeaturePlot(cluster3, features = c('EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4', # arterial
#                                    'NR2F2','EPHB4', 'CDH5', 'ESM1', 'FLT4', 'SOX18', 'SMARCA4')) # venous))
# dev.copy(png, file = "output-dissertation-figures_240720/240722-cluster3-arterial-venous-endo-markers-800x1095.png", width = 800, height = 1095)
# dev.off()
# 
# FeaturePlot(cluster5, features = c('EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4', # arterial
#                                    'NR2F2','EPHB4', 'CDH5', 'ESM1', 'FLT4', 'SOX18', 'SMARCA4')) # venous))
# dev.copy(png, file = "output-dissertation-figures_240720/240722-cluster5-arterial-venous-endo-markers-800x1095.png", width = 800, height = 1095)
# dev.off()

cluster3_5 <- subset(integrated, idents = c('3','5'))
FeaturePlot(cluster3_5, features = c('EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4', # arterial
                                     'NR2F2','EPHB4', 'CDH5', 'ESM1', 'FLT4', 'SOX18', 'SMARCA4')) # venous))
dev.copy(png, file = "output-dissertation-figures_240720/240722-cluster5-arterial-venous-endo-markers-800x1095.png", width = 800, height = 1095)
dev.off()

levels(integrated)

DotPlot(integrated, features = c('EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4', # arterial
                                 'SMARCA4','KRT18','PCDH12','CDH5', 'NR2F2','EPHB4', 'ESM1', 'FLT4', 'SOX18'), group.by = "seurat_clusters") +
  ggtitle("CELL MARKER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Features") +
  ylab("clusters")
dev.copy(png, file = "output-dissertation-figures_240720/240722-intergrated-arterial-venous-endo-markers-dotplot-500X800.png", width = 500, height = 800)
dev.off()

FeaturePlot(integrated, features = c('EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4', 'IGFBP5','FABP3',# arterial
                                 'PCDH12','CDH5', 'NR2F2','EPHB4', 'ESM1', 'FLT4', 'SOX18')) 
dev.copy(png, file = "output-dissertation-figures_240720/240723-intergrated-arterial-venous-endo-markers-featureplot-800x900.png", width = 800, height = 900)
dev.off()

DotPlot(integrated, features = c('EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4', 'IGFBP5','FABP3',# arterial
                                       'PCDH12','CDH5', 'NR2F2','EPHB4', 'ESM1', 'FLT4', 'SOX18'), group.by = "seurat_clusters") +
  ggtitle("CELL MARKER EXPRESSION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  coord_flip() +
  xlab("Features") +
  ylab("clusters")
dev.copy(png, file = "output-dissertation-figures_240720/240723-intergrated-arterial-venous-endo-markers-dotplot-500X800.png", width = 500, height = 800)
dev.off()

# = = = = = = = = = = = = = = = = =  = = = = = 230723 compare ECs from laminar (3+5) and static (2) =  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # 



# ======================== # MONOCOLE3 # ======================== # 
# SET SEURAT OBJ TO cell_data_set object for monocle3

cluster6_7[["UMAP"]] <- cluster6_7[["umap"]]


cds <- reduce_dimension(cluster6_7, reduction_method = "UMAP")
cds <- as.cell_data_set(cluster6_7)
cds
View(cds@reduce_dim_aux)

# to get cell metadata
colData(cds)

# to get gene metadata (the initial output: DataFrame with 30713 rows and 0 columns)
fData(cds)
rownames(fData(cds))
fData(cds)$gene_short_name <- rownames(fData(cds))
fData(cds)

# get gene count (rows=gene_short_name; column: cell ids (ie 'AAACCCTTcccc-2_1'))
counts(cds)

# add partition information (basically superclustes)
# assign all the cells to one partition
# create a name list/ vector: names = cell_ids; values: partition value = 1
# for all the row names that is the cells in the cell data, repeat it the value of one for the number of rows
reacreate.partition <- c(rep(1, length(cds@colData@rownames))) 
# assign to the names of cell ids (ie 'AAACCCTTcccc-2_1')
names(reacreate.partition) <- cds@colData@rownames 
# convert this as factor
recreate.partition <- as.factor(reacreate.partition)
# use this name list and assign to the partition
cds@clusters$UMAP$partitions <- recreate.partition



# assign cluster information

list_cluster <- cluster6_7@meta.data$seurat_clusters
cds@clusters$UMAP$clusters <- list_cluster



# add UMAP cell embeddinb/ UMAP coordinates
# to locate where the cell embeddings are stored (reducedDimNames)
cds
cds@int_colData@listData$reducedDims$UMAP <- cluster6_7@reductions$UMAP@cell.embeddings


# plot trajectory anlysis
# view all cluster
cluster.before.traj <- plot_cells(cds,
           color_cells_by = 'seurat_clusters',
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) +
  theme(legend.position = 'right')
# output = No trajectory to plot. Has learn_graph() been called yet? (same as the tutorial don't worry)

cluster.names <- plot_cells(cds,
           color_cells_by = "seurat_clusters",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c( 'grey','brown', 'cyan', 'violet', 'tan', 'yellow', 'green', 'brown', 'magenta', 'black')) +
  theme(legend.position = 'right')
           
combined <- cluster.before.traj+cluster.names
combined
dev.copy(png, file = "240707-static-UMAP-cluster-annotated-trajectory-2000x800.png", width = 2000, height = 900)
dev.off()
           

# learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)
plot_cells(cds,
           color_cells_by = 'seurat_clusters', label_cell_groups = F,
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5) + NoAxes() + NoLegend()
dev.copy(png, file = "output-dissertation-figures_240720//figure2/240806-cluster6_7-trajectory-288x332.png", width = 288, height = 332)
dev.off()


# order cells in pseudotime, you can pick the origin from here
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds@colData@listData$seurat_clusters == '0'))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE, 
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 7)
dev.copy(png, file = "240707-static-pseudotime-trajectory-1095x700.png", width = 1095, height = 700)
dev.off()

# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, 
       aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime, median), 
           fill = seurat_clusters)) +
  geom_boxplot()+theme_grey(base_size = 15)
dev.copy(png, file = "240707-static-annotated-pseudotime-trajectory-boxplot-1095x700.png", width = 1095, height = 700)
dev.off()

saveRDS(cluster6_7, file = '240707-static-only-for-trajectory-pilot-run.RDS')

# finding genes that chagne as a function of pseudotime, this process takes 1:20:00
deg <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg %>%
  arrage(q_value) %>%
  filter(status == 'OK') %>%
  head()





DotPlot(integrated, features = c(
                                 'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
                                 'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                                 'CD44','RUNX1','CD82','MYB', # EHT 
                                 'TAL1', 'PF4', # MEGA
                                 'AZU1', 'MPO', 'RNASE2', # GRAN
                                 'CD14','PDK4','NCF1', # MONOCYTE
                                 'C1QA','MRC1', # MACROPHAGE
                                 'HLA-DRA','HLA-DRB1'), # IMMUNO-CELLS
        group.by = 'Condition' )+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('')+
  ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/figure2/240812-intergrated-expression-distribution-dotplot-Con-1000X300.png", width = 1000, height = 300)
dev.off()


# COMPARE ARTERIAL ENDOTHLIUM IN EC CLUSTERS BY CONDITIONS
# CLUSTER3+5 <- LAMINAR
# CLUSTER 1 <- STATIC

Idents(integrated) <- 'seurat_clusters'
EC <- subset(integrated, ident = c('1','3','5'))

DimPlot(EC, group.by = 'seurat_clusters')
Idents(EC) <- 'Condition'
Idents(EC)

EC_con <- subset(EC, ident = c('Static','Laminar'))
DimPlot(EC_con, group.by = 'seurat_clusters')

#  DIP2A functions as a novel receptor that mediates the cardiovascularprotective effects of FSTL1.

DotPlot(EC_con, features = c('DIP2A','FSLT1'), # EC actuvation
        group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')

FeaturePlot(EC_con, feature = 'DIP2A')
# mechanosensitive markers
DotPlot(EC_con, features = c('KLF2','KLF4','NOS3', # REGULATE VASCULAR TONE
                             'CYP1A1','CYP1B1','PLPP3', # ANTI-ATHEROGENIC
                             'SLC9A3R2','PODXL','ADAM15', # VASCULAR HOMEOSTASIS AND EC SURVIVAL
                             'HMOX1','NQO1'), # EC actuvation
        group.by = 'Condition')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/figure2/240808-mechanosensitive-EC-800X300.png", width = 500, height =300)
dev.off()

FeaturePlot(integrated, features = c('KLF2','KLF4','NOS3', # REGULATE VASCULAR TONE
                                 'CYP1A1','CYP1B1','PLPP3', # ANTI-ATHEROGENIC
                                 'SLC9A3R2','PODXL','ADAM15', # VASCULAR HOMEOSTASIS AND EC SURVIVAL
                                 'HMOX1','NQO1'))

# laminar: #619CFF
# static: #00BA38
# FLOW-GENES:
# UP: 'KLF2','KLF4','NOS3', 'CYP1A1','CYP1B1','PLPP3','SLC9A3R2','PODXL','ADAM15', 'HMOX1','NQO1'
# DOWN: 'EDN1','APLN','ANGPT2','CITED2','DDAH1','THBS1','CAV1','SOD2'   ---> DDIN'T DOWNREGULATED EXCEPT EDN1
cluster_colour <- c('1' = '#D89000','3'='#39B600', '5' = '#39B600')

flowG <- c('PODXL','NQO1','SLC9A3R2','KLF2' ,'NOS3', 'ADAM15')
flow <- c('PODXL','NQO1','SLC9A3R2','KLF2' ,'NOS3', 'KLF4')
# REGULATE VASCULAR TONE and # VASCULAR HOMEOSTASIS AND EC SURVIVAL
VlnPlot(integrated, features = c(notch), ncol = 8,
        group.by = 'seurat_clusters',pt.size = 0) + scale_fill_manual(values = cluster_colour)
dev.copy(png, file = "output-dissertation-figures_240720/figure2/240815-KLF4.png", width = 960, height = 220)
dev.off()


# VlnPlot(EC_con, features = c('EDN1','APLN','ANGPT2','CITED2','DDAH1','THBS1','CAV1','SOD2'), 
#         group.by = 'Condition',
#         add.noise = F, pt.size = 0) # lol they aren't down-regulated

DimPlot(integrated)

# arterial marker --> no diff
VlnPlot(EC_con, features = c('NOTCH1','CXCR4','HEY1','GJA1','HES1', 'DLL4'), 
        group.by = 'seurat_clusters',
        add.noise = F, pt.size = 0)

VlnPlot(EC_con, features = c('CXCR4','DLL4','NOTCH1','NOTCH4','HEY1','HES1','GJA1'), pt.size= 0)
DotPlot(EC_con, features = c('CXCR4','DLL4','NOTCH1','NOTCH4','HEY1','HES1','GJA1'), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')

# CHECK NOTCH SINCE IT IS MECHANOSENSITIVE: it's activated in 1, 
Idents(EC_con) <- temp_ident
DotPlot(EC_con, features = c('DLL4','NOTCH1','NOTCH4','HEY1','HES1','GJA1'), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')


saveRDS(EC_con, flie = "RDS/240808-EC-static(1)-vs-laminar(35).RDS")

# to identify static EC
# PROLIFERATING: 'MKI67','TOP2A','BIRC5' (smh they are highly expressed in laminar but down in static)
# ARTERIAL: 'NOTCH1','ACKR3','CXCR4','HEY1','GJA1','HES1'
DimPlot(EC_con, label  =T, group.by = 'seurat_clusters')
DotPlot(EC_con, features = c('NOTCH1','ACKR3','CXCR4','HEY1','GJA1','HES1'))+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')

# notch: 'DLL4','NOTCH1','NOTCH4','HEY1','HES1','GJA1' --> all express
# don't use dot plot for this since it's the relative expression, use vlnplot instead
DotPlot(EC_con, features = c('DLL4','NOTCH1','NOTCH4','HEY1','HES1','GJA1'), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')

VlnPlot(EC_con,features = c('DLL4','NOTCH1','NOTCH4','HEY1','HES1','GJA1', # notch
                            'PGF','VEGFC','FLT1','KDR','FLT4','NRP1','NRP2', # vegf
                            'MKI67','TOP2A','BIRC5' ), #vegf
        group.by = 'seurat_clusters',
        pt.size = 0)

# VEGF: 'PGF','VEGFC','FLT1','KDR','FLT4','NRP1','NRP2'
VlnPlot(EC_con,features = c('HOPX', # ENDOTHEIUM
                            'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                            'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2'), #vegf
        group.by = 'seurat_clusters',
        pt.size = 0)

DotPlot(EC_con, features = c('PGF','VEGFC','FLT1','KDR','FLT4','NRP1','NRP2'))+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')

# plot findmarker
DotPlot(EC_con, features = c('DHRS3','SOX6','CXCR4','ACSM3','FAT3','ADGRG6','MMRN1','NLGN1',
                             'SLC9A3R2','TMEM173','A2M','KLF2','KLF4','FABP3','IGFBP5','APLNR'),

        group.by = 'Condition' )+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')
 dev.copy(png, file = "output-dissertation-figures_240720/figure2/240806-intergrated-expression-distribution-dotplot-Con-800x200.png", width = 800, height = 200)
dev.off()



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


# OXIDATIVE PHOSPHORYLATION + glycolysis
DotPlot(object = EC_con, features = c('ATP5MC1', 'COX5A', 'CYC1', #oxsphos
                                      'PDHA1', 'PKM', 'LDHA', 'HIF1A','TPI1')) + 
  ggtitle("OXIDATIVE PHOSPHORYLATION") +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(2, 7)) +
  RotatedAxis()+
  xlab("") +
  ylab("")

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

# # glyco-oxs
# # 'NRF2','MGST1','SOD2','ATP5MC1', 'COX5A', 'CYC1', 'PDHA1', 'PKM', 'LDHA', 'HIF1A','TPI1'
# # 
# # # FATTY ACID METABOLIC
# # 'SDHC', 'CA2', 'ADSL', 'ECHS1', 'DLD', 'ECI2','ODC1', 'IDH1', 'S100A10', 'CCDC58'
# # #  GLUTAMINE METABOLIC
# # 'GAD1', 'FAH', 'PYCR1', 'NIT2', 'GOT2', 'ADSL', 'GPT2', 'GLUD1' 
# 
# 
# DotPlot(EC_con, features = c('NRF2','MGST1','SOD2','COX5A', 'CYC1', 'PDHA1', 'PKM', 'LDHA', 'TPI1','HIF1A',
#                              'SDHC', 'CA2',  'DLD', 'ECHS1', 'ADSL','ODC1','ECI2', 'IDH1', 'S100A10', 'CCDC58',
#                              'GAD1', 'FAH', 'PYCR1', 'NIT2', 'GOT2', 'GPT2', 'GLUD1' ),
#         group.by = 'Condition' )+
#   ggtitle(NULL) +
#   scale_color_gradientn(colors = c('Blue', 'Red')) +
#   scale_size(range = c(1, 10)) +
#   RotatedAxis()+
#   xlab('') + ylab('')
# 
# 
# 
# 
# from the GSEA-result, get -NES == laminar-enriched
glycolysis_genes <- c('PGM1','GPI','PFKL','TPI1','GAPDH','PGK1', 'ENO1') # 8


tca_genes <- c('PDHA1','PDHB','DLAT','AKR1A1','DIP2A','ACSS2','ALDH2') # 8
  
# FROM SHINYGO OFF THE FINDMARKER RESULT FROM LAMINAR
oxidative_phosphorylation_genes <- c('SDHD','SDHB','COX5A','COX5B','COX7B','COX7C','COX17')

DotPlot(EC_con, features =c(glycolysis_genes, tca_genes, oxidative_phosphorylation_genes))+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/figure2/240809-ec_CON-glycolysis-TCA-OxPhos.png", width = 900, height = 300)
dev.off()

# glutathione metabolic process - 
#   "CHAC2", "GLRX2", "SOD1", "ETHE1", "IDH1", "GSTT2B", 
# "DPEP1", "SLC1A1", "G6PD", "GCLM", "SLC7A11", "GGT5"
# glua_genes <- c(
#   "CHAC2", "GLRX2", "SOD1", "ETHE1","SLC1A1", "G6PD", "GCLM", "SLC7A11") #8

# "FITM2", "HACD2", "ACSL4", "ELOVL6", "DGAT1", 
# "HSD17B12", "SLC25A1", "TECR", "HACD1", "ACOT7", 
# "NUDT8", "ELOVL1"

# lipid_metabolismgenes <- c(
#   "FITM2", "HACD2", "ACSL4", "ELOVL6", "DGAT1", 
#    "SLC25A1", "NUDT8", "ELOVL1")
glua_genes <- c(
 "GLRX2", "SOD1", "ETHE1", "G6PD", "GCLM", "SLC7A11")

lipid_metabolismgenes <- c(
  "HACD2", "ACSL4", "ELOVL6", 
  "SLC25A1", "NUDT8", "ELOVL1")

DotPlot(EC_con, features =c(glua_genes, lipid_metabolismgenes), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/figure6/240812-EC_CON-metabolism.png", width = 600, height = 260)
dev.off()

# glyco oxphos
DotPlot(EC_con, features =c(glycolysis_genes, tca_genes, oxidative_phosphorylation_genes), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/figure6/240812-oxphs.png", width = 900, height = 300)
dev.off()

# LAMINAR COMBINED
DotPlot(EC_con, features =c(glycolysis_genes, tca_genes, oxidative_phosphorylation_genes, glua_genes, lipid_metabolismgenes), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/figure6/240812-clusters135-all-remove4.png", width = 1000, height = 300)
dev.off()

# TGFb + WNT + HIPPO (YAP1)
# gsea-list checked

tgfb_pathway_genes <- c('ACVR1','SMAD2','SMAD3','SMAD4', 'ROCK1')
  # "SMAD2","SMAD3", "SMAD4", # TFs
  # 'PIK3CA','PIK3CB','PIK3R2', 'AKT1')# PI3K-AKT) 

# from GSEA ++ WNT specific TFs
wnt <- c('LRP6',"TCF7L1",  "SOX4", "FOXO1", "ZEB2", "GATA3")


# no on the GSEA list but the paper
# SAV1 is upstream of LAT1/2, activated LAT1/2 will stop YAP/TAZ from translocation (hippo on)
hippo_genes <- c("YAP1", "SAV1", "LATS1", "LATS2","TEAD1","TEAD2","TEAD3", "TEAD4")


DotPlot(EC_con, features =c(tgfb_pathway_genes, wnt, hippo_genes), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/240812-TGFB+WNT+hippo-EC-135clusters-820X262.png", width = 820, height = 262)
dev.off()

DotPlot(integrated, features =c(glycolysis_genes, tca_genes, oxidative_phosphorylation_genes, glua_genes, lipid_metabolismgenes), group.by = 'Condition')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('') + ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/figure2/240809-EC_CON-LAMINAR-ALL.png", width = 1400, height = 300)
dev.off()


# EC_con endo identity
DotPlot(EC_con, features = c(
                                 'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO
                                 'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
                                 'CD44','RUNX1','CD82','MYB', # EHT 
                                 'TAL1', 'PF4', # MEGA
                                 'AZU1', 'MPO', 'RNASE2', # GRAN
                                 'CD14','PDK4','NCF1', # MONOCYTE
                                 'C1QA','MRC1', # MACROPHAGE
                                 'HLA-DRA','HLA-DRB1'), # IMMUNO-CELLS
        group.by = 'Condition' )+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('')
ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/figure2/240806-intergrated-expression-distribution-dotplot-Con-1000X300.png", width = 1000, height = 300)
dev.off()

FeaturePlot(integrated, features = c('EFNB2', 'NOTCH1', 'DLL4', 'GJA5','HEY1', 'HEY2', 'SOX17', 'CXCR4', # arterial
                                     'NR2F2','EPHB4', 'CDH5', 'ESM1', 'FLT4', 'SOX18', 'SMARCA4')) # venous))

notch_r <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "JAG1", "JAG2", "DLL1", "DLL4")
notch_TF <- c("HES1", "HES5", "HES7", "HEY1", "HEY2", "HEYL", "RBPJ", "MAML1", "MAML2")

arterial_endo <- c("ESM1", "VEGFA", "ICAM1", "VCAM1", "NOS3", "THBD", "KDR", "MMP9", "HIF1A", "GATA2", "ECSCR", "JAG1")
Idents(EC_con)

DotPlot(EC_con, features = c('SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2',arterial_endo))+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('')
ylab('')


# arterial-specialised genes, NR2F2 is venous so should be low
# arterial_NOTCH <- c('ETS1', 'CXCR4','EFNB2','SOX7','SOX17','SOX18','DLL4','NOTCH1','NOTCH4','HEY1','HES1', 'CD93')
venous <- c('APLNR','FTH1','NRP2', 'KRT18','NR2F2','SMARCA4')

arterial <- c('GJA4','SOX7','SOX17','SOX18','DLL4', 'CD93')
notch <- c('JAG1','DLL4','NOTCH1','NOTCH4','HEY1','HEY2','HES1', 'GJA1','EFNB2')
# vegf <- c('PGF','VEGFC','FLT1','DKR','NRP1','NRP2')

DotPlot(EC_con, features = c(arterial, venous), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('')+ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/240815-arterial-venous-updated-660-300.png", width = 500, height = 300)
dev.off()

DotPlot(EC_con, features = c(notch), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('')+ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/240816-notch-updated-660-300.png", width = 600, height = 250)
dev.off()

# 'SOX17','HEY1','DLL4','GJA4','GJA5','EFNB2', # ARTERIAL ENDOTHELIUM
# 'KRT18','SMARCA4','NR2F2','PCDH12','FLT4','EPHB4', # VENOUS ENDO

FeaturePlot(EC_con, features = c(arterial, venous))
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/arterial-venous-dotplot-Con-1000X300.png", width = 2000, height = 2000)
dev.off()

DotPlot(integrated, features = c(arterial), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('')+ylab('')

DotPlot(EC_con, features = c(arterial,'GJA4','GJA5', venous,'KRT18','SMARCA4','PCDH12','FLT4'), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('')+ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/240812-arterial-venous-135clusters+HES1-820X262.png", width = 820, height = 262)
dev.off()

RidgePlot(EC_con, features = c(arterial,'GJA4','GJA5', venous,'KRT18','SMARCA4','PCDH12','FLT4'), group.by = 'seurat_clusters',ncol = 8)
VlnPlot(EC_con, features = c(arterial, venous), group.by = 'seurat_clusters',ncol = 8, pt.size = 0)


DotPlot(EC_con, features = c(glycolysis_genes, tca_genes, oxidative_phosphorylation_genes,arterial,'GJA4','GJA5', venous,'KRT18','SMARCA4','PCDH12','FLT4'), group.by = 'seurat_clusters')+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  xlab('')+ylab('')
dev.copy(png, file = "output-dissertation-figures_240720/EC comparison (static_vs_laminar)/240812-arterial-venous-135clusters+HES1-820X262.png", width = 820, height = 262)
dev.off()


DimPlot(EC_con, split.by = 'Condition', group.by = 'seurat_clusters') + NoAxes()

p <- DimPlot(EC_con, split.by = 'Condition', group.by = 'seurat_clusters') + NoAxes()

# Define the colors for specific clusters
p + scale_color_manual(values = c(
  '1' = '#F8766D' ,  # Cluster 1 to orange
  '3' = '#39B600'   ,# Cluster 3 to green
  '5' = '#39B600'    # Cluster 5 to green
))

# dotplot in EHT with ETS1 + notch to check if cluster EC primed towards EHT are arterialised HE
temp_idents <- cluster6_7@active.ident
cluster6_7@meta.data$seurat_clusters <- factor(cluster6_7@active.ident, levels = c('4','0','3','1','2'))
cluster6_7@meta.data$seurat_clusters <- factor(cluster6_7@active.ident, levels = c('2','1','3','0','4'))

DimPlot(cluster6_7, group.by = 'seurat_clusters')
DotPlot(EC_con, features = c('ETS1',notch),
        group.by = 'seurat_clusters' )+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  ylab(" ") +xlab("")

# ====
DotPlot(cluster6_7, features = c('ETS1',notch),
        group.by = 'seurat_clusters' )+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  ylab(" ") +xlab("")

# response to laminar fluid shear stress GSEA
NFE2L2/MAPK7/KLF2/KLF4
shear <- c('NFE2L2','MAPK7','KLF2','KLF4')
DotPlot(EC_con, features = c(shear),
        group.by = 'seurat_clusters' )+
  ggtitle(NULL) +
  scale_color_gradientn(colors = c('Blue', 'Red')) +
  scale_size(range = c(1, 10)) +
  RotatedAxis()+
  ylab(" ") +xlab("")

FeatureScatter(EC_con, 
               "JAG1", "DLL4", 
               group.by = "seurat_clusters", 
               pt.size = 3)

cc <- c('CDK4','CDK6',
        'CCNE2','CDK2',
        'CCND3',
        'CCNA2','CDK1',
        'CCNB1')

proliferation <- c('MKI67','CCNB2')
# CDKN1c: inhibit G1S transition
DotPlot(integrated, features = c('CDT1', cc, 'CDKN1C', proliferation, 'CD34'), group.by = 'Condition')+
       ggtitle(NULL) +
       scale_color_gradientn(colors = c('Blue', 'Red')) +
      scale_size(range = c(1, 10)) +
      RotatedAxis()+
      xlab('')+ylab('')
