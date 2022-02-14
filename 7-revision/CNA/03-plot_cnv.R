# This script plots the results obtained with inferCNV
# (https://github.com/broadinstitute/inferCNV/)


# Load packages
print("Loading packages...")
library(infercnv)
library(here)
library(glue)
set.seed(1234)


# Parse command line arguments
print("Parsing command-line arguments...")
args <- commandArgs(trailingOnly = TRUE)
donor_id <- args[[1]]


# Define paths
path_to_dir <- here()
path_to_obj <- glue("{path_to_dir}/results/inferCNV/infercnv_obj_{donor_id}.rds")
path_to_save <-  glue("{path_to_dir}/results/inferCNV/{donor_id}")
# path_to_save <- "~/Desktop/diogenes/"


# Read data
infercnv_obj <- readRDS(path_to_obj)


# Plot and save
plot_cnv(
  infercnv_obj = infercnv_obj,
  out_dir = path_to_save,
  title = "",
  obs_title = "",
  ref_title = "",
  cluster_references = FALSE,
  plot_chr_scale = TRUE,
  k_obs_groups = 1,
  custom_color_pal = color.palette(c("#f24153", "white", "#536db6"), between = c(2, 2)),
  color_safe_pal = TRUE,
  output_filename = glue("infercnv_custom_heatmap_{donor_id}"),
  output_format = "pdf",
  dynamic_resize = 0.25,
  write_expr_matrix = FALSE,
  useRaster = TRUE,
)
