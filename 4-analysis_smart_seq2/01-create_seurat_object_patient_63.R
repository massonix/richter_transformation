# This script creates a Seurat object specific for all cells from patient 063


# Load packages
library(Seurat)
library(tidyverse)


# Define paths
path_to_obj <- here::here("data/CLL_52/CLL_52.dgecounts.rds")
path_to_metadata <- here::here("data/CLL_52/CLL_52_metadata.csv")
path_to_patient_metadata <- here::here("data/sample_metadata.csv")
path_to_gene_names <- here::here("data/CLL_52/CLL_52.gene_names.txt")
path_to_63 <- here::here("results/R_objects/patient_63")
path_to_save <- str_c(path_to_63, "1.seurat_object_unfiltered.rds", sep = "/")


# Load data
counts <- readRDS(path_to_obj)
counts <- counts$readcount$exon$all
metadata <- read_csv(path_to_metadata)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$zUMI_barcode
gene_names_df <- as.data.frame(read_tsv(path_to_gene_names))
rownames(gene_names_df) <- gene_names_df$gene_id
gene_names_df <- gene_names_df[rownames(counts), ]
if (all(rownames(counts) == rownames(gene_names_df))) {
  rownames(counts) <- gene_names_df$gene_name
}
all(rownames(counts) %in% gene_names_df$gene_id)
metadata <- metadata[colnames(counts), ]
patient_metadata <- read_csv(path_to_patient_metadata)
patient_metadata <- patient_metadata[patient_metadata$donor_id == "63", ]
patient_metadata <- as.data.frame(patient_metadata)


# Create Seurat object
metadata_all <- left_join(metadata, patient_metadata, by = "sample_id")
metadata_all <- as.data.frame(metadata_all)
rownames(metadata_all) <- metadata_all$zUMI_barcode
if (all(rownames(metadata_all) == colnames(counts))) {
  seurat <- CreateSeuratObject(counts = counts, meta.data = metadata)
}


# Save
dir.create(path_to_63, showWarnings = FALSE)
saveRDS(seurat, path_to_save)