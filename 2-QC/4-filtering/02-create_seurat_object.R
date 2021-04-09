# This script loads the demultiplexed seurat object, the non-hashed expression
# matrices, and creates a single Seurat object


# Load packages
library(Seurat)
library(tidyverse)


# Define paths
path_to_hashed <- here::here("results/R_objects/1.seurat_hashed_demultiplexed.rds")
path_to_not_hashed <- here::here("1-cellranger_mapping/projects/BCLLATLAS_29/jobs/")
path_to_project_metadata <- here::here("1-cellranger_mapping/data/richter_metadata.csv")
path_to_scrublet <- here::here("results/tables/scrublet")
path_to_save <- here::here("results/R_objects/2.seurat_unfiltered.rds")


# Read data
metadata <- read_csv(path_to_project_metadata)
seurat_hashed <- readRDS(path_to_hashed)
files_to_load <- str_subset(
  list.dirs(path_to_not_hashed),
  "filtered_feature_bc_matrix"
)
gem_ids_not_hashed <- list.dirs(
  path = path_to_not_hashed,
  full.names = FALSE,
  recursive = FALSE
)
matrices_not_hashed <- purrr::map(files_to_load, Read10X)
seurat_list <- purrr::map2(matrices_not_hashed, gem_ids_not_hashed, function(mat, x) {
  seurat_obj <- CreateSeuratObject(counts = mat)
  seurat_obj$gem_id <- x
  new_barcodes <- str_c(x, colnames(seurat_obj), sep = "_")
  seurat_obj <- RenameCells(seurat_obj, new.names = new_barcodes)
  seurat_obj
})
names(seurat_list) <- gem_ids_not_hashed


# Add scrublet info to non-hashed cells
paths_to_scrublet <- list.files(path_to_scrublet, full.names = TRUE)
seurat_list <- purrr::map(gem_ids_not_hashed, function(x) {
  path <- str_subset(paths_to_scrublet, x)
  scrublet_df <- read_csv(path)
  scrublet_df$barcodes <- str_c(x, scrublet_df$barcodes, sep = "_")
  seurat_obj <- seurat_list[[x]]
  if (all(scrublet_df$barcodes == colnames(seurat_obj))) {
    warning("barcodes are equal")
    seurat_obj$scrublet_doublet_scores <- scrublet_df$scrublet_doublet_scores
    seurat_obj$scrublet_doublet_scores_scaled <- scale(
      scrublet_df$scrublet_doublet_scores,
      center = TRUE,
      scale = TRUE
    )
    seurat_obj$scrublet_predicted_doublet <- scrublet_df$scrublet_predicted_doublet
    return(seurat_obj)
    
  } else{
    warning("barcodes are not equal")
    return(NULL)
  }
})
names(seurat_list) <- gem_ids_not_hashed


# Merge not hashed objects
seurat_not_hashed <- seurat_list[[1]]
for (i in 2:length(seurat_list)) {
  print(gem_ids_not_hashed[i])
  seurat_not_hashed <- merge(x = seurat_not_hashed, y = seurat_list[[i]])
  seurat_list[[i]] <- NA
}
rm(seurat_list)


# Remove assay containing hashtag oligonucleotide expression (HTO)
seurat_hashed[["HTO"]] <- NULL


# Set hashing-related variables to NA
seurat_hashed$nCount_HTO <- NULL
seurat_hashed$nFeature_HTO <- NULL
seurat_hashed$orig.ident <- NULL
seurat_not_hashed$orig.ident <- NULL
seurat_not_hashed$HTO_maxID <- "NA"
seurat_not_hashed$HTO_secondID <- "NA"
seurat_not_hashed$HTO_margin <- "NA"
seurat_not_hashed$HTO_classification <- "NA"
seurat_not_hashed$HTO_classification.global <- "NA"
seurat_not_hashed$hash.ID <- "NA"
seurat_not_hashed$hashing_snr <- NA


# Include donor id and library name
metadata <- filter(metadata, type != "hashed_hto")
library_names <- metadata$library_name
donor_ids <- metadata$donor_id
names(library_names) <- metadata$gem_id
names(donor_ids) <- metadata$gem_id
seurat_hashed$library_name <- library_names[seurat_hashed$gem_id]
seurat_hashed$donor_id <- donor_ids[seurat_hashed$gem_id]
seurat_not_hashed$library_name <- library_names[seurat_not_hashed$gem_id]
seurat_not_hashed$donor_id <- donor_ids[seurat_not_hashed$gem_id]


# Merge
seurat_hashed$is_hashed <- TRUE
seurat_not_hashed$is_hashed <- FALSE
seurat_not_hashed@meta.data <- seurat_not_hashed@meta.data[, colnames(seurat_hashed@meta.data)]
seurat <- merge(x = seurat_hashed, y = seurat_not_hashed)


# Save
saveRDS(seurat, path_to_save)
