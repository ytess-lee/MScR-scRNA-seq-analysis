# MiloR
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(tidyverse)
library(patchwork)
library(monocle3)
library(SeuratObject)
library(glmGamPoi)
library(miloR)

# set project directory
setwd("F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis")

integrated <- readRDS("RDS/240719-integrated-with-cc-data.RDS") 

# example 
data(sim_trajectory)
milo.meta <- sim_trajectory$meta
milo.obj <- Milo(sim_trajectory$SCE)
milo.obj
milo.obj <- buildGraph(milo.obj, k=20, d=30)
milo.obj <- makeNhoods(milo.obj, k=20, d=30, refined=TRUE, prop=0.2)
milo.obj <- calcNhoodDistance(milo.obj, d=30)
milo.obj <- countCells(milo.obj, samples="Sample", meta.data=milo.meta)

milo.design <- as.data.frame(xtabs(~ Condition + Sample, data=milo.meta))
milo.design <- milo.design[milo.design$Freq > 0, ]
rownames(milo.design) <- milo.design$Sample
milo.design <- milo.design[colnames(nhoodCounts(milo.obj)),]

milo.res <- testNhoods(milo.obj, design=~Condition, design.df=milo.design)
head(milo.res)

# subset the seurat obj by conditions to reduce processing time
static <- subset(integrated, idents = 'Static')
prestim <- subset(integrated, idents = 'Prestimulation')

# convert seurat obj to singlecellexperimetn (SCE) obj
sce <- as.SingleCellExperiment(static)

# create milo obj from sce obj
milo_static <- Milo(sce)

# build a KNN graph and count neighorhoods
# Build kNN graph (using precomputed neighbors)
milo_static <- buildGraph(milo_static, k = 20, d = 30)  # Adjust k and d based on your data

# Count cells in neighborhoods
milo_static <- countCells(sample =milo_static)
