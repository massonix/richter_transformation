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
