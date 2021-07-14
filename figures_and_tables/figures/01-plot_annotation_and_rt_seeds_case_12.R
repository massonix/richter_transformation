# This script plots the UMAPs of the annotation and RT seed cells for the 5
# patients


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)


# Source scripts and variables
source(here::here("bin/utils.R"))


# Load Seurat object
path_to_12 <- here::here("results/R_objects/6.seurat_annotated_12.rds")
seurat <- readRDS(path_to_12)


# Individual panels
umap_annotation <- plot_annotation(
  seurat_obj = seurat,
  pt_size = 0.5,
  colors_reference = color_annotations,
  patient_id = "12",
  nothing = FALSE
)
umap_annotation <- umap_annotation +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    legend.position = c(0.05, 0.8),
    legend.text = element_text(size = 9)
  )

genes_interest <- c("CXCR4", "CD24", "CD27", "MIR155HG", "CCND2", "TOP2A",
                    "PCNA", "MZB1", "IGHM", "XBP1")
dot_plot <- plot_dot_plot(
  seurat,
  goi = genes_interest,
  colors_reference = color_annotations,
  patient_id = "12"
)
dot_plot <- dot_plot +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
violin_plot_s_phase <- plot_violin_plot(
  seurat_obj = seurat,
  continuous_var = "S.Score",
  ylab = "S Phase Score",
  colors_reference = color_annotations,
  patient_id = "12"
)
violin_plot_g2m_phase <- plot_violin_plot(
  seurat_obj = seurat,
  continuous_var = "G2M.Score",
  ylab = "G2M Phase Score",
  colors_reference = color_annotations,
  patient_id = "12"
)