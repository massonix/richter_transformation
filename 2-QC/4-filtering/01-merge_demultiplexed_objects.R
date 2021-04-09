# This script creates a single Seurat object containing cells from all samples


# Load packages
library(Seurat)
library(tidyverse)


# Define paths
path_to_demultiplexed <- here::here("results/R_objects/demultiplexed/")
path_to_save <- here::here("results/R_objects/1.seurat_hashed_demultiplexed.rds")


# Read data
files_to_load <- list.files(path_to_demultiplexed)
gem_ids <- files_to_load %>%
  str_remove("seurat_") %>%
  str_remove("_demultiplexed.rds")
files_to_load2 <- str_c(path_to_demultiplexed, files_to_load, sep = "")
seurat_list <- purrr::map2(files_to_load2, gem_ids, function(path, x) {
  seurat_obj <- readRDS(path)
  seurat_obj$gem_id <- x
  seurat_obj
})
names(seurat_list) <- gem_ids


# Merge Seurat objects
print("Merging Seurat objects...")
seurat_hashed <- seurat_list[[1]]
for (i in 2:length(seurat_list)) {
  print(gem_ids[i])
  seurat_hashed <- merge(x = tonsil_hashed, y = seurat_list[[i]])
  seurat_list[[i]] <- NA
}
rm(seurat_list)


# Save
print("Saving...")
saveRDS(seurat_hashed, path_to_save)
