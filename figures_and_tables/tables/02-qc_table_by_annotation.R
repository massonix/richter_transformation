# This script generates a quality control (QC) table for each time-point and
# cluster


# Load packages
library(tidyverse)
library(Seurat)


# Define parameters
path_12 <- here::here("results/R_objects/6.seurat_annotated_12.rds")
path_19 <- here::here("results/R_objects/6.seurat_annotated_19.rds")
path_63 <- here::here("results/R_objects/patient_63/3.seurat_annotated.rds")
path_365 <- here::here("results/R_objects/6.seurat_annotated_365.rds")
path_3299 <- here::here("results/R_objects/6.seurat_annotated_3299.rds")


# Load Seurat objects
seurat_12 <- readRDS(path_12)
seurat_19 <- readRDS(path_19)
seurat_63 <- readRDS(path_63)
seurat_365 <- readRDS(path_365)
seurat_3299 <- readRDS(path_3299)


# Create list
seurat_list <- list(
  "12" = seurat_12,
  "19" = seurat_19,
  "63" = seurat_63,
  "365" = seurat_365,
  "3299" = seurat_3299
)
rm(seurat_12, seurat_19, seurat_365, seurat_3299, seurat_63)


# Make table
n_cells_tables <- purrr::map(seurat_list, function(seurat_obj) {
  df <- seurat_obj@meta.data %>%
    group_by(time_point, sample_description_FN, subproject, factor(annotation_final)) %>%
    dplyr::count(.drop = FALSE, )
  df
})
n_cells_table <- bind_rows(n_cells_tables, .id = "case")


# Save
colnames(n_cells_table) <- c(
  "Case",
  "Time point",
  "Description",
  "Subproject",
  "Annotation",
  "Num of cells passing filters"
)
# Save
openxlsx::write.xlsx(
  n_cells_table,
  "results/tables/paper/scRNA-seq_qc_table_by_annotation.xlsx"
)
write_delim(
  n_cells_table,
  file = "results/tables/paper/scRNA-seq_qc_table_by_annotation.csv",
  delim = ";"
)
