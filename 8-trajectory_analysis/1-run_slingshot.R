# Load packages
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)


# Load data
path_to_obj <- here::here("results/R_objects/patient_63/7.seurat_list_annotated_reprocessed.rds")
seurat_list <- readRDS(path_to_obj)


# Subset to CLL clusters only
seurat_list$`12`$annotation_final <- Idents(seurat_list$`12`)
seurat_list$`19`$annotation_final <- Idents(seurat_list$`19`)
seurat_list$`3299`$annotation_final <- Idents(seurat_list$`3299`)
seurat_list$`12` <- subset(
  seurat_list$`12`,
  idents = c("CXCR4-CD27+", "CXCR4+CD27-", "MIR155HG+")
)
seurat_list$`19` <- subset(
  seurat_list$`19`,
  idents = c("CXCR4-CD27+", "CXCR4+CD27-", "MIR155HGhi CLL-like")
)
seurat_list$`3299` <- subset(
  seurat_list$`3299`,
  idents = c("CD27+S100A4+", "CD83hiMIR155HGhi", "CXCR4-CD27+", "CXCR4+CD27-", "CD83loMIR155HGhi")
)


# Merge
donors <- names(seurat_list)
seurat <- merge(
  x = seurat_list$`12`,
  y = c(seurat_list$`19`, seurat_list$`3299`)
)
seurat <- seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2500) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(dims = 1:30, reduction = "pca", group.by.vars = "donor_id") %>%
  RunUMAP(reduction = "harmony", dims = 1:30)


DimPlot(seurat, group.by = "donor_id", pt.size = 0.75)
FeaturePlot(seurat, "CXCR4", pt.size = 0.5) +
  scale_color_viridis_c(option = "inferno")
FeaturePlot(seurat, "CD27", pt.size = 0.5) +
  scale_color_viridis_c(option = "inferno")
FeaturePlot(seurat, "MIR155HG", pt.size = 0.5) +
  scale_color_viridis_c(option = "inferno")
FeaturePlot(seurat, "PCDH9", pt.size = 0.5) +
  scale_color_viridis_c(option = "inferno")



# Run Slingshot
seurat <- FindNeighbors(seurat, dims = 1:30, reduction = "harmony")
seurat <- FindClusters(seurat, resolution = 0.6)
DimPlot(seurat)
seurat <- ScaleData(seurat, features = rownames(seurat))
sce <- as.SingleCellExperiment(seurat)

n_cells <- Matrix::rowSums(seurat[["RNA"]]@counts > 0)
gene_qc <- n_cells %>% 
  as.data.frame() %>% 
  ggplot(aes(n_cells)) + 
  geom_histogram(bins = 100, alpha = 0.75) +
  scale_x_log10("Number of cells") +
  theme_bw() 
gene_qc
sce <- sce[names(n_cells)[n_cells > 30], ]
sce <- slingshot(
  sce,
  clusterLabels = "seurat_clusters",
  reducedDim =  "HARMONY"
)
seurat$slingshot_pseudotime1 <- sce$slingPseudotime_1
seurat$slingshot_pseudotime2 <- sce$slingPseudotime_2
FeaturePlot(seurat, features = "slingshot_pseudotime2", pt.size = 0.35) +
  scale_color_viridis_c(option = "magma")


# Run trade-Seq
# http://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html
seurat2 <- FindVariableFeatures(seurat, nfeatures = 4000, selection.method = "dispersion")
sce <- tradeSeq::fitGAM(sce, genes = VariableFeatures(seurat2), verbose = TRUE)
ATres <- associationTest(sce)


# Save
saveRDS(sce, "results/R_objects/SingleCellExperiment_slingshot.rds")
saveRDS(seurat, "results/R_objects/Seurat_slingshot.rds")






###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Do the same for patient 12
seurat_list$`12`$annotation_final <- Idents(seurat_list$`12`)
seurat <- seurat_list$`12`
seurat <- seurat %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = rownames(seurat)) %>%
  RunPCA(features = VariableFeatures(seurat)) %>%
  RunUMAP(dims = 1:20, reduction = "pca")
seurat <- subset(
  seurat,
  idents = c("CXCR4-CD27+", "CXCR4+CD27-", "MIR155HG+")
)

# Filter out genes
n_cells <- Matrix::rowSums(seurat[["RNA"]]@counts > 0)
gene_qc <- n_cells %>% 
  as.data.frame() %>% 
  ggplot(aes(n_cells)) + 
  geom_histogram(bins = 100, alpha = 0.75) +
  scale_x_log10("Number of cells") +
  theme_bw() 
gene_qc

sce <- as.SingleCellExperiment(seurat)
sce <- sce[names(n_cells)[n_cells > 15], ]
sce <- slingshot(
  sce,
  clusterLabels = "annotation_final",
  reducedDim =  "PCA"
)
seurat$slingshot_pseudotime <- sce$slingPseudotime_1
FeaturePlot(seurat, features = "slingshot_pseudotime", pt.size = 1) +
  scale_color_viridis_c(option = "magma")


# TradeSeq
seurat2 <- FindVariableFeatures(seurat, nfeatures = 6000)
input_genes <- VariableFeatures(seurat2)[VariableFeatures(seurat2) %in% rownames(sce)]
sce <- tradeSeq::fitGAM(
  sce,
  genes = input_genes,
  verbose = TRUE
)
ATres <- associationTest(sce)


oxphos_genes <- AnnotationDbi::select(
  x = org.Hs.eg.db,
  keys = "GO:0006119",
  keytype = "GO",
  columns = "SYMBOL"
)$SYMBOL
seurat <- AddModuleScore(seurat,features = list(oxphos_genes), name = "oxphos")
FeaturePlot(seurat, features = "oxphos1") +
  scale_color_viridis_c(option = "magma")


topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])
