# This script tries to model the LN-PB recirculation of CLL cells using Monocle3


# Load packages
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(monocle3)
library(tidyverse)
library(plotly)


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
umap_df <- as.data.frame(seurat@reductions$umap@cell.embeddings)
umap_df$annotation_final <- Idents(seurat)
umap_gg <- plot_ly(
  umap_df,
  x = ~UMAP_1,
  y = ~UMAP_2,
  z = ~UMAP_3,
  color = ~factor(annotation_final),
  # colors = color_palette,
  size = 0.001,
  alpha = 0.4
)
umap_gg <- umap_gg %>%
  add_markers()
umap_gg
DimPlot(seurat)


# Run Monocle
cds <- as.cell_data_set(seurat)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(
  cds,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
# Choose the root cell based on the highest expression of CXCR4
root <- CellSelector(FeaturePlot(seurat, features = "CXCR4", reduction = "pca", dims = 1:2))
cds <- order_cells(cds, root_cells = root)
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
seurat2 <- as.Seurat(cds, assay = "RNA")
seurat$monocle3_pseudotime <- seurat2$monocle3_pseudotime
umap_df$pseudotime <- seurat$monocle3_pseudotime
umap_pseudotime_gg <- plot_ly(
  umap_df,
  x = ~UMAP_1,
  y = ~UMAP_2,
  z = ~UMAP_3,
  color = ~pseudotime,
  # colors = color_palette,
  size = 0.001,
  alpha = 0.4
)
umap_pseudotime_gg <- umap_pseudotime_gg %>%
  add_markers()
umap_pseudotime_gg



# # Plot expression of genes as a function of pseudotime
# selected_genes <- c("CXCR4", "MIR155HG", "CD27")
# plots <- purrr::map(selected_genes, function(x) {
#   p <- seurat@meta.data %>%
#     mutate(gene = seurat[["RNA"]]@data[x, ]) %>% 
#     ggplot(aes_string("monocle3_pseudotime", "gene")) +
#       geom_point() +
#       theme_bw()
#   p
# })
# plots

################################################################################
################################################################################

seurat <- readRDS("results/R_objects/Seurat_slingshot.rds")
DimPlot(seurat)
cds <- as.cell_data_set(seurat)
cds <- align_cds(cds, alignment_group = "donor_id")
cds <- reduce_dimension(cds)
cds <- learn_graph(cds = cds, use_partition = TRUE)
plot_cells(cds, color_cells_by = "seurat_clusters")

####################################################

seurat <- readRDS("results/R_objects/Seurat_slingshot.rds")
seurat_list <- SplitObject(seurat, split.by = "donor_id")
for (i in seq_along(seurat_list)) {
  seurat_list[[i]] <- seurat_list[[i]] %>%
    NormalizeData() %>%
    FindVariableFeatures()
}
features <- SelectIntegrationFeatures(seurat_list)
for (i in seq_along(along.with = seurat_list)) {
  seurat_list[[i]] <- seurat_list[[i]] %>%
    ScaleData(features = features) %>%
    RunPCA(features = features)
}
anchors <- FindIntegrationAnchors(seurat_list, reduction = "rpca", dims = 1:30)
integrated <- IntegrateData(anchors, dims = 1:30)



integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated)
DimPlot(integrated)


cds <- as.cell_data_set(integrated)
cds <- cluster_cells(cds, resolution = 0.5)
cds <- learn_graph(cds)
plot_cells(
  cds,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)
