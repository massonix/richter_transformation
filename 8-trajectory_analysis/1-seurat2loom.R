# This script tries to model the LN-PB recirculation of CLL cells using Monocle3


# Load packages
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(slingshot)
library(tidyverse)


# Load data
path_to_obj <- here::here("results/R_objects/patient_63/7.seurat_list_annotated_reprocessed.rds")
seurat_list <- readRDS(path_to_obj)
seurat <- seurat_list$`12`
seurat$annotation_final <- Idents(seurat)
rm(seurat_list)


# Source script with functions
source(here::here("bin/utils.R"))


# Keep CLL cells only and reprocess
selected_clusters <- c("CXCR4-CD27+", "CXCR4+CD27-", "MIR155HG+")
seurat <- subset(seurat, idents = selected_clusters)


# Reprocess Seurat object
seurat <- seurat %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "time_point", dims = 1:25, reduction = "pca") %>%
  RunUMAP(reduction = "harmony", dim = 1:25, n.components = 3)
seurat$UMAP_1 <- seurat@reductions$umap@cell.embeddings[, "UMAP_1"]
seurat$UMAP_2 <- seurat@reductions$umap@cell.embeddings[, "UMAP_2"]
seurat$UMAP_3 <- seurat@reductions$umap@cell.embeddings[, "UMAP_3"]


# Cluster
seurat <- FindNeighbors(seurat, reduction = "harmony", dims = 1:25)
seurat<- FindClusters(seurat, resolution = 0.5)
DimPlot(seurat)
# seurat <- FindSubCluster(
#   seurat,
#   cluster = "5",
#   graph.name = "RNA_snn",
#   resolution = 0.2,
#   subcluster.name = "mir155"
# )
# DimPlot(seurat, group.by = "mir155")
# seurat$clusters_input_slingshot <- seurat$mir155
# seurat$clusters_input_slingshot <- factor(seurat$clusters_input_slingshot)
# levels(seurat$clusters_input_slingshot) <- as.character(
#   0:length(levels(seurat$clusters_input_slingshot))
# )
# DimPlot(seurat, group.by = "clusters_input_slingshot")


# Run slingshot
sce <- as.SingleCellExperiment(seurat)
sce <- slingshot(
  sce,
  clusterLabels = "seurat_clusters",
  reducedDim =  "HARMONY"
)
seurat$slingshot_pseudotime <- sce$slingPseudotime_1
FeaturePlot(seurat, features = "slingshot_pseudotime", pt.size = 1) +
  scale_color_viridis_c(option = "viridis")


# Run tradeSeq

# Save
saveRDS(seurat, here::here("results/R_objects/seurat_12_trajectory.rds"))
