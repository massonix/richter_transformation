# Again, CreateChromatinAssay needs to download something from the Internet,
# which cannot be done in the computing nodes of the cluster. Thus, we will
# isolate this part of the script.


# Load packages
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(GEOquery)
library(tidyverse)


# Define paths
path_to_annotation <- here::here("7-revision/Penter2021/scATAC-seq/tmp/genomic_annotation.rds")
path_to_mat <- here::here("7-revision/Penter2021/scATAC-seq/tmp/mat_tmp.rds")
path_to_peaks <- here::here("7-revision/Penter2021/scATAC-seq/tmp/combined_peaks_tmp.rds")
path_to_frags <- here::here("7-revision/Penter2021/scATAC-seq/tmp/frag_objs_tmp.rds")
path_to_gse <- here::here("results/R_objects/gse_atac_obj.rds")
path_to_save <- here::here("results/R_objects/Penter2021/scATAC-seq/1.seurat_obj_init_Penter2021_atac.rds")


# Read data
annotation <- readRDS(path_to_annotation)
mat <- readRDS(path_to_mat)
combined_peaks <- readRDS(path_to_peaks)
frag_objs <- readRDS(path_to_frags)


# Initialize Seurat object
chrom_assay <- CreateChromatinAssay(
  counts = mat,
  ranges = combined_peaks,
  genome = "hg38",
  fragments = frag_objs,
  annotation = annotation,
  sep = c("-", "-"),
  min.features = 200
)
seurat_obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")


# Add metadata
seurat_obj$dataset_id <- str_sub(colnames(seurat_obj), start = 1, end = 10)
seurat_obj$cell_barcode <- colnames(seurat_obj)
gse <- readRDS(path_to_gse)
gsm_list <- GSMList(gse)
gsm_list <- gsm_list[unique(seurat_obj$dataset_id)]
metadata_list <- purrr::map(gsm_list, function(x) {
  gsm <- x@header
  df <- data.frame(
    dataset_id = gsm$geo_accession,
    library_strategy = gsm$library_strategy,
    sample_description = gsm$title
  )
  df$tissue <- gsm$characteristics_ch1 %>%
    str_subset("tissue") %>%
    str_remove("tissue: ")
  df$wbc <- gsm$characteristics_ch1 %>%
    str_subset("wbc") %>%
    str_remove("wbc: ") %>%
    as.numeric()
  df$disease_state <- gsm$characteristics_ch1 %>%
    str_subset("disease state") %>%
    str_remove("disease state: ")
  tumor_purity <- gsm$characteristics_ch1 %>%
    str_subset("cd19+") %>%
    str_split(":")
  df$tumor_purity <- as.numeric(str_remove(tumor_purity[[1]][2], "^ "))
  df
})
metadata_df <- bind_rows(metadata_list)
new_metadata <- left_join(seurat_obj@meta.data, metadata_df, by = "dataset_id")
rownames(new_metadata) <- new_metadata$cell_barcode
seurat_obj@meta.data <- new_metadata


# Save
saveRDS(seurat_obj, path_to_save)
