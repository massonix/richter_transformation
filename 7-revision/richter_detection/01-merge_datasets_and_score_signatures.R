# This script integrates the training and test sets to detect
# Richter cells in unseen peripheral blood samples
# Thus, we will only work with peripheral blood samples


# Load packages
print("Loading packages...")
library(Seurat)
library(harmony)
library(tidyverse)
library(readxl)
library(here)
library(UCell)


# Define parameters
print("Defining parameters...")
path_to_training_set <- here::here("results/R_objects/Penter2021/3.seurat_obj_clustered_Penter2021.rds")
path_to_test_set <- here::here("results/R_objects/6.seurat_annotated_12.rds")
path_to_cll9 <- here::here("results/R_objects/Penter2021/4.seurat_obj_CLL9_Penter2021.rds")
path_to_save <- here::here("results/R_objects/seurat_training_and_test.rds")
path_to_gsea <- here::here("results/tables/Penter2021/gsea_RT_vs_CLL_Penter2021.xlsx")


# Source functions
print("Sourcing functions...")
source(here::here("bin/utils.R"))


# Read data
print("Reading data...")
training_set <- readRDS(path_to_training_set)
test_set <- readRDS(path_to_test_set)
cll9 <- readRDS(path_to_cll9)
gsea <- read_excel(path_to_gsea)


# Subset to keep only peripheral blood
print("Subsetting...")
training_set <- subset(training_set, tissue == "Peripheral blood")
training_set$donor_id <- str_extract(training_set$sample_id, "CLL.")
Idents(training_set) <- "donor_id"
training_set <- subset(training_set, downsample = 5000)
training_set


# Homogenize metadata and define labels
print("Homogenizing metadata...")
training_set$is_rt <- ifelse(
  colnames(training_set) %in% colnames(cll9)[cll9$annotation == "RT"],
  "RT",
  "CLL"
)
test_set$is_rt <- ifelse(
  str_detect(test_set$annotation_final, "RT"),
  "RT",
  "CLL"
)
training_set$is_training <- "training"
test_set$is_training <- "test"
merged <- merge(x = training_set, y = test_set)


# Dimensionality reduction and integration
print("Performing dimensionality reduction...")
hvg <- find_assay_specific_features(
  merged,
  assay_var = "is_training",
  n_features = 5000
)
merged <- integrate_assays(
  merged,
  assay_specific = TRUE,
  shared_hvg = hvg,
  assay_var = "is_training",
  n_dim = 30
)
merged <- RunUMAP(merged, dims = 1:30, reduction = "harmony")


# Calculate OXPHOS and BCR scores
print("Calculating scores...")
oxphos_related <- c(
  "ATP metabolic process",                                        
  "oxidative phosphorylation",                                   
  "aerobic respiration",
  "mitochondrial translation",
  "NADH dehydrogenase complex assembly",
  "mitochondrial transmembrane transport",
  "cellular detoxification"
)
bcr_related <- c(
  "regulation of B cell activation",
  "response to calcium ion",
  "humoral immune response",
  "complement activation, classical pathway",
  "defense response to bacterium",
  "phagocytosis, recognition",
  "adaptive immune response",
  "humoral immune response mediated by circulating immunoglobulin",
  "B cell receptor signaling pathway",
  "phagocytosis, engulfment",
  "plasma membrane invagination"
)
gsea <- as.data.frame(gsea)
oxphos_signature <- gsea[gsea$Description %in% oxphos_related, "core_enrichment", drop = TRUE]
oxphos_signature <- str_split(oxphos_signature, pattern = "/")
oxphos_signature <- unique(unlist(oxphos_signature))
bcr_signature <- gsea[gsea$Description %in% bcr_related, "core_enrichment", drop = TRUE]
bcr_signature <- str_split(bcr_signature, pattern = "/")
bcr_signature <- unique(unlist(bcr_signature))
merged <- AddModuleScore_UCell(
  merged,
  features = list(
    oxphos_score = oxphos_signature,
    bcr_score = bcr_signature
  )
)


# Save
print("Saving...")
saveRDS(merged, path_to_save)

