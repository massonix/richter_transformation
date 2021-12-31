# This script creates the accessibility matrix and Seurat object for
# Penter et al. 2021


# Load packages
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(GEOquery)
library(tidyverse)
library(glue)


# Define paths and parameters
frag_path <- "~/Desktop/richter_transformation/data/scATAC/scATAC-seq"
macs_path <- "/Users/rmassoni/opt/miniconda3/bin/macs2"
path_to_gse <- "/Users/rmassoni/Desktop/richter_transformation/data/scATAC/gse_atac_obj.rds"
# path_to_gse <- here::here("results/R_objects/gse_atac_obj.rds")
n_barcodes <- 30000
# n_barcodes <- 900

# Call cells
gsm_dirs <- frag_path %>%
  list.dirs(full.names = FALSE) %>%
  str_subset("^GSM")
n_frag_dfs <- purrr::map(gsm_dirs, function(x) {
  path <- list.files(
    glue("{frag_path}/{x}"),
    full.names = TRUE,
    pattern = "fragments_per_cell.txt"
  )
  df <- read_delim(path, col_names = c("cell_barcode", "n_fragments"))
  df
})
names(n_frag_dfs) <- gsm_dirs
n_frag_df <- bind_rows(n_frag_dfs, .id = "dataset_id")
n_frag_df <- n_frag_df %>%
  arrange(desc(n_fragments)) %>%
  mutate(position = 1:nrow(n_frag_df))
(n_frag_gg <- n_frag_df %>%
  ggplot(aes(position, n_fragments)) +
    geom_line() +
    labs(x = "barcodes", y = "Number of Fragments") +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    geom_vline(xintercept = n_barcodes, linetype = "dashed", color = "red"))
sel_cells <- n_frag_df %>%
  filter(position <= n_barcodes) %>%
  select("dataset_id", "cell_barcode")


# Call peaks
frag_objs <- purrr::map(gsm_dirs, function(x) {
  path <- list.files(
    glue("{frag_path}/{x}"),
    full.names = TRUE,
    pattern = "_fragments.tsv.gz$"
  )
  cells <- sel_cells[sel_cells$dataset_id == x, "cell_barcode", drop = TRUE]
  frag_obj <- CreateFragmentObject(path, cells = cells)
  frag_obj
})
names(frag_objs) <- gsm_dirs
peaks_list <- purrr::map(frag_objs, function(x) {
  peaks <- CallPeaks(x, macs2.path = macs_path)
  peaks
})
combined_peaks <- GenomicRanges::reduce(c(
  peaks_list$GSM4982189,
  peaks_list$GSM4982190,
  peaks_list$GSM4982191,
  peaks_list$GSM4982192,
  peaks_list$GSM4982193,
  peaks_list$GSM4982194,
  peaks_list$GSM4982195,
  peaks_list$GSM4982196,
  peaks_list$GSM4982197
))
combined_peaks <- keepStandardChromosomes(
  combined_peaks,
  pruning.mode = "coarse"
)
combined_peaks <- subsetByOverlaps(
  x = combined_peaks,
  ranges = blacklist_hg38_unified,
  invert = TRUE
)
peak_widths <- width(combined_peaks)
combined_peaks <- combined_peaks[peak_widths  < 10000 & peak_widths > 20]


# Create accessibility matrix
mat <- FeatureMatrix(
  fragments = frag_objs,
  features = combined_peaks,
  cells = sel_cells$cell_barcode
)


# Set up Seurat object
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
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
# saveRDS(seurat_obj, "~/Desktop/richter_transformation/results/R_objects/seurat_atac.rds")

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




# See https://github.com/timoast/signac/issues/156
# Corre la següent línia:
# gzip -dc GSM4982189_CLL9_1.fragments.tsv.gz | awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,"GSM4982189_"$4,$5}'