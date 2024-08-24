# Milo implementation using graph based features for neighbourhood sampling and SpatialFDR
# https://github.com/MarioniLab/miloR/issues/293
library(ggplot2)
library(Seurat)
# library(SeuratWrappers)
library(dplyr)
library(tidyverse)
library(patchwork)
library(SeuratObject)
library(miloR)
library(SingleCellExperiment)
library(singleCellTK)
Sys.setenv(LANG = "en")
# ## Load dummy data
# data(sim_trajectory)
# milo.meta <- sim_trajectory$meta
# milo.obj <- Milo(sim_trajectory$SCE)
# 
# ## Build KNN graph neighbourhoods
# milo.obj <- buildGraph(milo.obj, k=40, d=13)
# milo.obj <- makeNhoods(milo.obj, k=40, d=13, refined=TRUE, prop=0.1, refinement_scheme="graph")
# 
# ## Count cells in nhoods
# milo.obj <- countCells(milo.obj, samples="Sample", meta.data=milo.meta)
# 
# ## Test for differential abundance
# milo.design <- as.data.frame(xtabs(~ Condition + Sample, data=milo.meta))
# milo.design <- milo.design[milo.design$Freq > 0, ]
# rownames(milo.design) <- milo.design$Sample
# milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]
# milo_res <- testNhoods(milo.obj, design=~Condition, design.df=milo.design, fdr.weighting="graph-overlap")

# set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")

integrated <- readRDS("RDS/240719-integrated-with-cc-data.RDS") 
Idents(integrated) <- "Condition" 

# load milo obj
# Convert from Seurat object to SCE object 
integrated.sce <- as.SingleCellExperiment(integrated)

# create milo obj from sce
milo.obj <- Milo(integrated.sce)

# build knn 
milo.obj <- buildGraph(milo.obj, k=40, d=13)
milo.obj <- makeNhoods(milo.obj, k=40, d=13, refined=TRUE, prop=0.1, refinement_scheme="graph")

# check distribution
plotNhoodSizeHist(milo.obj)

# create df
Meta.data <- data.frame(colData(milo.obj))
# # Create a new column 'sample' by concatenating 'ident' and 'dataset' columns
Meta.data$sample <- paste(Meta.data$ident, Meta.data$dataset, sep = "_")

# counting cells in neighbourhoods
milo.obj <-countCells(milo.obj, meta.data = Meta.data, samples="sample")

# test for DA
in.design <- Meta.data[,c("sample", "Condition")]
in.design <- distinct(in.design)
rownames(in.design) <- in.design$sample

## Test for differential abundance
# reorder rownames to match columns of nhoodCounts(milo)
in.design <- in.design[colnames(nhoodCounts(milo.obj)), , drop = FALSE]
in.design

# save RDS
saveRDS(in.design, file = 'F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/milo-rds/240803-milo-script-in_desgin.RDS')
saveRDS(milo.obj, file = 'F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/milo-rds/240803-milo-script-integrated_milo.RDS')

# define our experimental design
rownames(in.design) <- in.design$sample
da_results <- testNhoods(milo.obj, design = ~ Condition, design.df = in.design, fdr.weighting = "graph-overlap")

da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

saveRDS(da_results, file = 'F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/milo-rds/240803-milo-script(speedo)-DA_results.RDS')


# visualisation
milo.obj <- readRDS(file = 'F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/milo-rds/240803-milo-script-integrated_milo.RDS')
da_results <- readRDS(file = 'F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/milo-rds/240803-milo-script(speedo)-DA_results.RDS')

write.csv(da_results, file = 'F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/milo-rds/da_results.csv')

milo.obj <- buildNhoodGraph(milo.obj)
plotUMAP(milo.obj)+
  plot_layout(guides="collect")
dev.copy(png, file = "240804-miloR(speeded)-milo-obj(Nhood)-collect-2000x2000.png", width = 2000, height = 2000)
dev.off()

plotNhoodGraphDA(milo.obj, da_results, alpha=0.05, res_column = "logFC") +
  plot_layout(guides="collect")
dev.copy(png, file = "240804-miloR(speeded)-milo-obj+DA_result-res-added-1400x1000.png", width = 1400, height = 1000)
dev.off()

# 240815 - FDR =5% --> white
# it's weird it wouldn't plot against logFC even if i specified, and csv looks fine i wonder if it's a computational power issue
# okay i tried wtih pvalue, same issue
# annotateNhoods(x, da.res, coldata_col, subset.nhoods = NULL)
annotated_da_res <- annotateNhoods(milo.obj, da_results, coldata_col='seurat_clusters')
write.csv(annotated_da_res, "240814-annotated_da_res.csv", row.names = FALSE)


annotated_da_res <- read.csv(file = "240814-annotated_da_res.csv")

annotated_da_res$seurat_clusters <- factor(annotated_da_res$seurat_clusters, levels = 0:9)

plotDAbeeswarm(annotated_da_res, group.by = 'seurat_clusters', alpha = 0.1, subset.nhoods = NULL) + 
  scale_color_manual(values = c("#8cc2e4", "#90c47d", "#af8fd2", "#e48484", "#e4ad84", "#e4e484", "#8ce1e1", "#e48398", "#c2c2c2", "#deb698", 'white'))

# plotDAbeeswarm(annotated_da_res, group.by = 'seurat_clusters', alpha = 0.1, subset.nhoods = NULL)+
#   scale_colour_gradient(low = "blue", high = "red")
# 
# plotDAbeeswarm(annotated_da_res, group.by = "seurat_clusters") +
#   scale_fill_gradient2(low = "red", mid = "white",
#                        high = "blue", midpoint = 0, 
#                        na.value = "grey50", guide = "colourbar")


# Normal beeswarm plot using ggbeeswarm -> this is madness
ggplot(annotated_da_res, aes(x = logFC, y = seurat_clusters , color = seurat_clusters)) +
  geom_beeswarm(alpha = 0.1, size = 1.5) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal() +
  labs(x = "LogFC", y = "seurat_clusters", title = "Beeswarm Plot of Log Fold Change by Seurat Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# dotplot
ggplot(annotated_da_res, aes(x =logFC , y = factor(seurat_clusters), color = logFC)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(x = "Log Fold Change", y = "Seurat Clusters", title = "Dot Plot of Log Fold Change by Seurat Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




milo.obj <- buildNhoodGraph(milo.obj)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo.obj, dimred = "umap", colour_by="Condition", text_by = "seurat_clusters", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(embryo_milo, da_results, layout="umap",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")









# calculate enrichment manually -> nope

# Define the cluster-to-condition mapping
prestimulation_clusters <- c(0, 1, 4)
laminar_clusters <- c(3, 5)
static_clusters <- c(2, 6, 7, 8, 9)

# Assign conditions based on seurat_clusters
annotated_da_res$condition <- ifelse(annotated_da_res$seurat_clusters %in% prestimulation_clusters, "Prestimulation",
                              ifelse(annotated_da_res$seurat_clusters %in% laminar_clusters, "Laminar",
                                     ifelse(annotated_da_res$seurat_clusters %in% static_clusters, "Static", NA)))


# Assuming cell_data is your dataframe with neighborhood and condition information -> it's bad
prop_table <- table(annotated_da_res$Nhood, annotated_da_res$condition)
props <- prop.table(prop_table, margin = 1)

# Example: Ratio of Static vs Prestimulation
enrichment_static_vs_prestim <- props[, "Static"] / props[, "Prestimulation"]

# Example: Ratio of Laminar vs Prestimulation
enrichment_laminar_vs_prestim <- props[, "Laminar"] / props[, "Prestimulation"]

# Add to your DA results
annotated_da_res$Enrichment_Static_vs_Prestim <- enrichment_static_vs_prestim[as.character(annotated_da_res$Nhood)]
annotated_da_res$Enrichment_Laminar_vs_Prestim <- enrichment_laminar_vs_prestim[as.character(annotated_da_res$Nhood)]
head(props)
