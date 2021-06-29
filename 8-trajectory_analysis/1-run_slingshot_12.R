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
seurat <- subset(seurat, slingshot_pseudotime < 35)
markers_dfs <- purrr::map(c("CXCR4", "CD27", "MIR155HG"), function(x) {
  df <- data.frame(
    pseudotime = seurat$slingshot_pseudotime,
    expression = seurat[["RNA"]]@scale.data[x, ],
    gene = x
  )
  df
})
markers_df <- bind_rows(markers_dfs)
markers_gg <- markers_df %>%
  ggplot(aes(pseudotime, expression, color = gene)) +
  # geom_point(size = 0.1) +
  geom_smooth() +
  theme_bw()
markers_gg


# Differential Expression Analysis
clusters <- c("CXCR4-CD27+", "CXCR4+CD27-", "MIR155HG+")
dea_list <- purrr::map(clusters, function(x) {
  df <- FindMarkers(
    seurat,
    ident.1 = x,
    only.pos = FALSE,
    logfc.threshold = 0.2
  )
  df <- df %>%
    rownames_to_column("gene") %>%
    filter(p_val_adj < 0.001) %>%
    arrange(desc(avg_log2FC))
})


# GSEA
set.seed(1234)
gsea <- purrr::map(dea_list, function(df) {
  gene_list <- df$avg_log2FC
  names(gene_list) <- df$gene
  gsea_results <- gseGO(
    gene_list,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    minGSSize = 25,
    maxGSSize = 300,
    seed = TRUE
  )
  gsea_results <- clusterProfiler::simplify(gsea_results, cutoff = 0.35)
  gsea_results@result <- gsea_results@result %>%
    arrange(desc(NES))
  gsea_results
})
names(gsea) <- clusters
DT::datatable(gsea$`CXCR4-CD27+`@result)
DT::datatable(gsea$`CXCR4+CD27-`@result)
DT::datatable(gsea$`MIR155HG+`@result)


signatures <- list(
  b_cell_signaling = gsea$`CXCR4-CD27+`@geneSets$`GO:0050853`,
  oxphos = gsea$`CXCR4-CD27+`@geneSets$`GO:0006119`,
  mitochondrial_translation = gsea$`CXCR4-CD27+`@geneSets$`GO:0032543`
)

seurat <- AddModuleScore(
  seurat,
  features = signatures,
  name = names(signatures)
)
seurat

signatures_pots <- purrr::map(
  c("b_cell_signaling1", "oxphos2", "mitochondrial_translation3"),
  function(x) {
    p1 <- FeaturePlot(seurat, x, pt.size = 0.7) +
      scale_color_viridis_c(option = "magma")
    p1
  })
signatures_pots


####






# TradeSeq
seurat2 <- FindVariableFeatures(seurat, nfeatures = 6000)
input_genes <- VariableFeatures(seurat2)[VariableFeatures(seurat2) %in% rownames(sce)]
sce <- tradeSeq::fitGAM(
  sce,
  genes = input_genes,
  verbose = TRUE
)
ATres <- associationTest(sce)
selected_rows <- 
  !is.na(ATres$waldStat) &
  !is.na(ATres$pvalue) &
  !is.na(ATres$meanLogFC)
ATres_sub <- ATres[selected_rows, ]
ATres_sub <- rownames_to_column(ATres_sub, "gene")
ATres_sub <- arrange(ATres_sub, pvalue, desc(waldStat))
ATres_sub$pvalue_adjusted <- p.adjust(ATres_sub$pvalue, method = "fdr")
DT::datatable(ATres_sub)



topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

