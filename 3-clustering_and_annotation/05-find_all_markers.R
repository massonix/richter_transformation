# This script finds and saves the markers for each of the clusters
# annotated in the previous notebook


# Load packages
library(Seurat)
library(tidyverse)


# Define paths
path_to_obj <- here::here("results/R_objects/patient_63/7.seurat_list_annotated_reprocessed.rds")
path_to_save <- here::here("results/tables/markers/")


# Load data
seurat_list <- readRDS(path_to_obj)


# Find Markers
markers_all <- purrr::map(seurat_list, function(seurat_obj) {
  df <- FindAllMarkers(
    seurat_obj,
    logfc.threshold = 0.65,
    only.pos = TRUE,
    test.use = "wilcox"
  )
  df
})
markers_all_dfs <- purrr::map(markers_all, function(df) {
  clusters <- unique(df$cluster)
  dfs <- purrr::map(clusters, function(x) {
    df_sub <- df[df$cluster == x & df$p_val_adj < 0.01, ]
    df_sub <- arrange(df_sub, desc(avg_log2FC))
    df_sub <- df_sub[, c(7, 1, 5, 2, 3, 4, 6)]
    df_sub
  })
  names(dfs) <- clusters
  dfs
})


# Gene Ontology Enrichment Analysis
enrichr_outputs_all <- purrr::map(markers_all_dfs, function(dfs) {
  enrichr_outputs <- purrr::map(dfs, function(df) {
    out <- enrichr(
      genes = df$gene,
      databases = "GO_Biological_Process_2018"
    )$GO_Biological_Process_2018
    selected_terms <- out$Adjusted.P.value < 0.01 & out$Odds.Ratio > 2.5
    out <- out[selected_terms, ]
    out
  })
  enrichr_outputs
})


# Save
for (patient in names(markers_all_dfs)) {
  path_to_save_markers <- str_c(
    path_to_save,
    "markers_annotated_clusters_patient_",
    patient,
    ".xlsx"
  )
  path_to_save_GO <- str_c(
    path_to_save,
    "GO_terms_annotated_clusters_patient_",
    patient,
    ".xlsx"
  )
  openxlsx::write.xlsx(markers_all_dfs[[patient]], file = path_to_save_markers)
  openxlsx::write.xlsx(enrichr_outputs_all[[patient]], file = path_to_save_GO)
}


