# This script creates a single excel file with the differential expresion analysis
# (RT vs CLL) and gene set enrichment analysis (GSEA) results for our reanalysis
# of the scRNA-seq dataset of Penter et al. (2021).


# Load packages
library(tidyverse)
library(Seurat)
library(openxlsx)
library(here)


# Define parameters
path_to_dea <- here("results/R_objects/Penter2021/dea_RT_vs_CLL_Penter2021.rds")
path_to_gsea <- here("results/R_objects/Penter2021/gsea_RT_vs_CLL_Penter2021_object.rds")
path_to_save <- here("results/tables/paper/DEA_and_GSEA_Penter_et_al.xlsx")
alpha <- 0.05

# Read data
dea <- readRDS(path_to_dea)
gsea <- readRDS(path_to_gsea)


# Subset GSEA
gsea_df <- gsea@result %>%
  dplyr::filter(p.adjust < alpha) %>%
  dplyr::arrange(desc(NES))


# Save
output_list <- list(
  RT_vs_CLL = dea,
  GSEA = gsea_df
)
openxlsx::write.xlsx(output_list, path_to_save)