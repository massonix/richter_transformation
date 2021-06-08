# This script tries to model the LN-PB recirculation of CLL cells using Monocle3


# Load packages
library(Seurat)
library(SeuratWrappers)
library(monocle3)
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
seurat <- process_seurat(seurat, dims = 1:25)
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
root <- CellSelector(DimPlot(seurat))
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


# Plot expression of genes as a function of pseudotime
selected_genes <- c("CXCR4", "MIR155HG", "CD27")
plots <- purrr::map(selected_genes, function(x) {
  p <- seurat@meta.data %>%
    mutate(gene = seurat[["RNA"]]@data[x, ]) %>% 
    ggplot(aes_string("monocle3_pseudotime", "gene")) +
      geom_point() +
      theme_bw()
  p
})
plots
