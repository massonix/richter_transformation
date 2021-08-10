# This script produces the two panels for the main figure describing the mutually
# exclusive pattern between OXPHOS and BCR


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)
library(ggridges)
library(ggforce)


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
umap_oxphos <- FeaturePlot(seurat, features = "oxphos1", pt.size = 0.05)
umap_oxphos <- umap_oxphos +
  ggtitle("OXPHOS Score") +
  scale_color_viridis_c(option = "magma") +
  theme(
    plot.title = element_text(face = "plain", size = 9),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    legend.position = c(0, 0.75),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.25, "cm")
  )
umap_bcr <- FeaturePlot(seurat, features = "bcr_signaling2", pt.size = 0.05)
umap_bcr <- umap_bcr +
  ggtitle("BCR Signaling Score") +
  scale_color_viridis_c(option = "magma") +
  theme(
    plot.title = element_text(face = "plain", size = 9),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.line = element_blank(),
    legend.position = c(0, 0.75),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm")
  )
umaps <- (umap_oxphos / umap_bcr)


# Ridge plots
ridgeplot_oxphos <- plot_ridge(
  seurat,
  feature = "oxphos1",
  x_lab = "OXPHOS Score",
  colors_reference = color_annotations,
  patient_id = "12"
)
ridgeplot_oxphos <- ridgeplot_oxphos +
  theme(
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 8)
  )
ridgeplot_bcr <- plot_ridge(
  seurat,
  feature = "bcr_signaling2",
  x_lab = "BCR Signaling Score",
  colors_reference = color_annotations,
  patient_id = "12"
)
ridgeplot_bcr <- ridgeplot_bcr +
  theme(
    axis.text = element_text(size = 5.5),
    axis.title = element_text(size = 8)
  )

# Violin plots
## Define RT seed cells
seurat$is_rt <- ifelse(str_detect(seurat$annotation_final, "RT"), "RT", "CLL")
seurat$sample_description_RM <- str_c(
  seurat$time_point,
  seurat$sample_description,
  sep = "_"
)
color_violin <- c("#dcdddc", "#cd8899")
violin_oxphos <- seurat@meta.data %>%
  ggplot(aes(sample_description_RM, oxphos1, fill = is_rt, color = is_rt)) +
    geom_violin(alpha = 0.35, position = "dodge") +
    ggforce::geom_sina(size = 0.5, alpha = 1) +
    scale_color_manual(values = color_violin) +
    scale_fill_manual(values = color_violin) +
    labs(x = "", y = "OXPHOS Score", color = "", fill = "") +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
violin_oxphos <- violin_oxphos + theme(axis.text.x = element_blank())
violin_bcr <- seurat@meta.data %>%
  ggplot(aes(sample_description_RM, bcr_signaling2, fill = is_rt, color = is_rt)) +
    geom_violin(alpha = 0.35) +
    ggforce::geom_sina(size = 0.5, alpha = 1) +
    scale_color_manual(values = color_violin) +
    scale_fill_manual(values = color_violin) +
    labs(x = "", y = "BCR Signaling Score", color = "", fill = "") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0.95, color = "black", size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      legend.text = element_text(size = 7)
    )


# Arrange figure
fig <- 
  (umap_oxphos / umap_bcr) |
  (ridgeplot_oxphos / ridgeplot_bcr) |
  (violin_oxphos/ violin_bcr)
fig


# Save
ggsave(
  filename = here::here("results/plots/paper/02-rt_oxphos_and_bcr_scRNA.pdf"),
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 11, 
  units = "cm"
)
