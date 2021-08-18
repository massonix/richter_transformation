# This script loads the patient-specific cluster markers and creates a single
# excel file (supplementary table)


# Load packages
library(tidyverse)
library(Seurat)


# Define parameters
path_12 <- here::here("results/tables/markers/markers_annotated_clusters_patient_12.rds")
path_19 <- here::here("results/tables/markers/markers_annotated_clusters_patient_19.rds")
path_63 <- here::here("results/tables/markers/markers_annotated_clusters_patient_63.rds")
path_365 <- here::here("results/tables/markers/markers_annotated_clusters_patient_365.rds")
path_3299 <- here::here("results/tables/markers/markers_annotated_clusters_patient_3299.rds")
paths_vec <- c("12" = path_12, "19" = path_19, "63" = path_63, "365" = path_365, "3299" = path_3299)


# Read data
markers_list <- purrr::map(paths_vec, readRDS)


# Save
openxlsx::write.xlsx(
  markers_list,
  here::here("results/tables/paper/scRNA-seq_markers_annotated_clusters.xlsx")
)
