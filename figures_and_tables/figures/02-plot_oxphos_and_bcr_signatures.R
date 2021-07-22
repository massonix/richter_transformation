# This script produces the two panels for the main figure describing the mutually
# exclusive pattern between OXPHOS and BCR


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)


# Source scripts and variables
source(here::here("bin/utils.R"))


# Load Seurat object and gene sets
path_to_12 <- here::here("results/R_objects/6.seurat_annotated_12.rds")
seurat <- readRDS(path_to_12)
gene_sets_list <- readRDS(here::here("6-differential_expression_analysis/tmp/gene_sets_gsea_analysis.rds"))


# Calculate signatures
signatures <- list(
  oxphos = gene_sets_list$`GO:0006119`,
  bcr_signaling = gene_sets_list$`GO:0050853`
)
mt_genes <- str_c("MT-", signatures$oxphos, sep = "") %in% rownames(seurat)
signatures$oxphos[mt_genes] <- str_c("MT-", signatures$oxphos[mt_genes], sep = "")
seurat <- AddModuleScore(seurat, features = signatures, name = names(signatures))


# Project signatures UMAP
umaps_signatures <- FeaturePlot(
  seurat,
  features = c("oxphos1", "bcr_signaling2"), 
  cols = c(),
  blend = TRUE,
  blend.threshold = 0.5
)
umaps_signatures[[1]]

