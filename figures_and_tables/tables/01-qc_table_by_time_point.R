# This script generates a quality control (QC) table for each time-point and
# patient.


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
qc_tables <- purrr::map(seurat_list, function(seurat_obj) {
  df <- seurat_obj@meta.data %>%
    group_by(time_point, sample_description_FN, subproject, tissue) %>%
    summarise(
      n_cells_passing_filters = n(),
      mean_library_size = round(mean(nCount_RNA), 3),
      mean_n_detected_genes = round(mean(nFeature_RNA), 3),
      mean_pct_mt = round(mean(pct_mt), 3)
    )
  df
})
qc_table <- bind_rows(qc_tables, .id = "case")
qc_table$is_hashed <- ifelse(
  qc_table$subproject == "BCLLATLAS_10",
  "TRUE",
  "FALSE"
)
qc_table$method <- ifelse(
  qc_table$subproject == "CLL_52",
  "Smart-seq2",
  "Chromium v3"
)


# Save
colnames(qc_table) <- c(
  "Case",
  "Time point",
  "Description",
  "Subproject",
  "Tissue",
  "Num of cells passing filters",
  "Mean library size (total UMI or counts)",
  "Mean num of detected genes",
  "Mean mitochondrial expression (%)",
  "Is hashed",
  "Method"
)
reordered_colnames <- c(
  "Case",
  "Time point",
  "Description",
  "Subproject",
  "Tissue",
  "Method",
  "Is hashed",
  "Num of cells passing filters",
  "Mean library size (total UMI or counts)",
  "Mean num of detected genes",
  "Mean mitochondrial expression (%)"
)
qc_table <- qc_table[, reordered_colnames]


# Save
openxlsx::write.xlsx(
  qc_table,
  "results/tables/paper/scRNA-seq_qc_table_by_time_point.xlsx"
)
write_delim(
  qc_table,
  file = "results/tables/paper/scRNA-seq_qc_table_by_time_point.csv",
  delim = ";"
)
