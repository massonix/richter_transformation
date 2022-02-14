# This script runs inverCNV (https://github.com/broadinstitute/inferCNV/)


# Load packages
print("Loading packages...")
library(Seurat)
library(here)
library(glue)
library(tidyverse)
library(infercnv)
set.seed(1234)


# Parse command line arguments
print("Parsing command-line arguments...")
args <- commandArgs(trailingOnly = TRUE)
donor_id <- args[[1]]


# Define parameters
# path_to_dir <- "~/Desktop/CNAG/mnt_clust/RICHTER/current"
print("Defining parameters...")
path_to_dir <- here::here()
path_to_utils <- glue::glue("{path_to_dir}/bin/utils.R")
path_to_obj <- ifelse(
  donor_id == "63",
  glue::glue("{path_to_dir}/results/R_objects/patient_63/3.seurat_annotated.rds"),
  glue::glue("{path_to_dir}/results/R_objects/6.seurat_annotated_{donor_id}.rds")
)
path_to_gene_order_f <- glue::glue("{path_to_dir}/7-revision/CNA/tmp/gencode_v21_gen_pos.complete.txt")
path_to_gene_order_donor <- glue::glue("{path_to_dir}/7-revision/CNA/tmp/gencode_v21_gen_pos.complete_{donor_id}.txt")
path_to_cell_annot <- glue::glue("{path_to_dir}/7-revision/CNA/tmp/cell_annotation_{donor_id}.txt")
path_to_outdir <- glue::glue("{path_to_dir}/results/inferCNV/donor_{donor_id}")
path_to_save <- glue::glue("{path_to_dir}/results/inferCNV/infercnv_obj_{donor_id}.rds")
dir.create(path_to_outdir, recursive = TRUE)

cutoff_infercnv <- if (donor_id == "63") 1 else 0.1


# Source functions
print("Sourcing functions...")
source(path_to_utils)


# Load data
print("Loading data...")
seurat <- readRDS(path_to_obj)
gene_order_file <- read_tsv(
  path_to_gene_order_f,
  col_names = c("gene", "chromosome", "start", "end")
)


# Write cell annotations file
print("Writing cell annotations file...")
seurat$cell_barcodes <- colnames(seurat)
# rt_vs_cll <- rt_vs_cll_dict[[donor_id]]
# seurat$annotation_infercnv <- rt_vs_cll[seurat$annotation_final]
seurat$annotation_infercnv <- ifelse(
  str_detect(seurat$annotation_final, "RT"),
  "malignant",
  "normal"
)
seurat$annotation_infercnv[seurat$annotation_infercnv == "malignant"] <- str_c(
  seurat$annotation_infercnv[seurat$annotation_infercnv == "malignant"],
  seurat$sample_description_FN[seurat$annotation_infercnv == "malignant"],
  sep = "_"
)
downsampled_cll <- sample(
  colnames(seurat)[seurat$annotation_infercnv == "normal"],
  sum(str_detect(seurat$annotation_infercnv, "malignant")),
  replace = FALSE
)
selected_barcodes <- c(
  colnames(seurat)[str_detect(seurat$annotation_infercnv, "malignant")],
  downsampled_cll
)
seurat <- subset(seurat, cells = selected_barcodes)
cell_annotations <- seurat@meta.data[, c("cell_barcodes", "annotation_infercnv")]
write_tsv(cell_annotations, path_to_cell_annot, col_names = FALSE)


# Obtain gene ordering file
print("Preparing and saving gene annotation file...")
gene_order_file <- as.data.frame(gene_order_file)
gene_order_file$gene <- str_remove(gene_order_file$gene, "\\|ENSG.*$")
seurat <- subset(
  seurat,
  features = rownames(seurat)[rownames(seurat) %in% gene_order_file$gene]
)
gene_order_file <- gene_order_file[match(rownames(seurat), gene_order_file$gene), ]
write_tsv(gene_order_file, path_to_gene_order_donor, col_names = FALSE)


# Create inferCNV object
print("Creating inferCNV object...")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = seurat[["RNA"]]@counts,
  annotations_file = path_to_cell_annot,
  delim = "\t",
  gene_order_file = path_to_gene_order_donor,
  ref_group_names = c("normal")
)


# Run inferCNV
print("Running inferCNV...")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = cutoff_infercnv,
  out_dir = path_to_outdir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE
)


# Save
print("Saving...")
saveRDS(infercnv_obj, path_to_save)
