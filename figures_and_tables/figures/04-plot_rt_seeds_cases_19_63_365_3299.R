# This script plots the UMAPs of early RT seeding cells for all
# cases except 12


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)


# Source scripts and variables
source(here::here("bin/utils.R"))


# Load Seurat object
path_to_19 <- here::here("results/R_objects/6.seurat_annotated_19.rds")
path_to_63 <- here::here("results/R_objects/patient_63/3.seurat_annotated.rds")
path_to_365 <- here::here("results/R_objects/6.seurat_annotated_365.rds")
path_to_3299 <- here::here("results/R_objects/6.seurat_annotated_3299.rds")
paths_to_load <- c(
  "19" = path_to_19,
  "63" = path_to_63,
  "365" = path_to_365,
  "3299" = path_to_3299
)
seurat_list <- purrr::map(paths_to_load, readRDS)


# UMAP RT seeds
umaps_seed_cells <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, x) {
  p <- plot_split_annotation(
    seurat_obj,
    pt_size = 0.8,
    split_by = "time_point",
    colors_reference = color_annotations,
    patient_id = x,
    n_col = 3
  )
  p <- p &
    NoLegend() &
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank()
    )
  p
})


fig <- (
  umaps_seed_cells$`19` /
    umaps_seed_cells$`63` /
    umaps_seed_cells$`365` /
    umaps_seed_cells$`3299` 
)
fig <- fig +
  plot_layout(heights = c(1.5, 1, 1, 1))


# Save
ggsave(
  filename = here::here("results/plots/paper/rt_seed_cells_supplementary.pdf"),
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 28, 
  units = "cm"
)



