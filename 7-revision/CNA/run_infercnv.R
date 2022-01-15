# This script runs inverCNV (https://github.com/broadinstitute/inferCNV/)


# Load packages
library(Seurat)
library(here)
library(glue)
library(tidyverse)
library(infercnv)
set.seed(1234)


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
donor_id <- args[[1]]


# Define parameters
# path_to_dir <- here::here()
path_to_dir <- "~/Desktop/CNAG/mnt_clust/RICHTER/current"
path_to_utils <- glue::glue("{path_to_dir}/bin/utils.R")
path_to_obj <- glue::glue("{path_to_dir}/results/R_objects/6.seurat_annotated_{donor_id}.rds")
path_to_gene_order_f <- glue::glue("{path_to_dir}/7-revision/CNA/tmp/gencode_v21_gen_pos.complete.txt")
path_to_cell_annot <- glue::glue("{path_to_dir}/7-revision/CNA/tmp/cell_annotation_{donor_id}.txt")
path_to_outdir <- glue::glue("{path_to_dir}/results/inferCNV/donor_{donor_id}")
path_to_save <- glue::glue("{path_to_dir}/results/inferCNV/infercnv_obj_{donor_id}.rds")
dir.create(path_to_outdir, recursive = TRUE)

cutoff_infercnv <- if (donor_id == "63") 1 else 0.1


# Source functions
source(path_to_utils)



# Load data
# seurat <- readRDS(here::here("results/R_objects/6.seurat_annotated_12.rds"))
seurat <- readRDS(path_to_obj)


# Create inferCNV object
seurat$cell_barcodes <- colnames(seurat)
rt_vs_cll <- rt_vs_cll_dict[[donor_id]]
seurat$annotation_infercnv <- rt_vs_cll[seurat$annotation_final]
downsampled_cll <- sample(
  colnames(seurat)[seurat$annotation_infercnv == "normal"],
  sum(seurat$annotation_infercnv == glue("malignant_{donor_id}")),
  replace = FALSE
)
selected_barcodes <- c(
  colnames(seurat)[seurat$annotation_infercnv == glue("malignant_{donor_id}")],
  downsampled_cll
)
seurat <- subset(seurat, cells = selected_barcodes)
cell_annotations <- seurat@meta.data[, c("cell_barcodes", "annotation_infercnv")]
write_tsv(cell_annotations, path_to_cell_annot, col_names = FALSE)
gene_order_file <- read_tsv(
  "https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v21_gen_pos.complete.txt",
  col_names = c("gene", "chromosome", "start", "end")
)
gene_order_file <- as.data.frame(gene_order_file)
gene_order_file$gene <- str_remove(gene_order_file$gene, "\\|ENSG.*$")
seurat <- subset(
  seurat,
  features = rownames(seurat)[rownames(seurat) %in% gene_order_file$gene]
)
gene_order_file <- gene_order_file[match(rownames(seurat), gene_order_file$gene), ]
write_tsv(gene_order_file, path_to_gene_order_f, col_names = FALSE)
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = seurat[["RNA"]]@counts,
  annotations_file = path_to_cell_annot,
  delim = "\t",
  gene_order_file = path_to_gene_order_f,
  ref_group_names = c("normal")
)


# Run inferCNV
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = cutoff_infercnv,
  out_dir = path_to_outdir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE
)


# SaveRDS
saveRDS(infercnv_obj, path_to_save)
