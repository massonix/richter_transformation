# This script plots the supplementary figure with the OXPHOS/BCR results


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)
library(ggrepel)
library(ggpubr)
library(gridGraphics)


# Source scripts and variables
source(here::here("bin/utils.R"))


# Load data
path_to_19 <- here::here("results/R_objects/6.seurat_annotated_19.rds")
path_to_63 <- here::here("results/R_objects/patient_63/3.seurat_annotated.rds")
path_to_365 <- here::here("results/R_objects/6.seurat_annotated_365.rds")
path_to_3299 <- here::here("results/R_objects/6.seurat_annotated_3299.rds")
paths_to_load <- c(
  "19" = path_to_19,
  "365" = path_to_365,
  "3299" = path_to_3299,
  "63" = path_to_63
)
seurat_list <- purrr::map(paths_to_load, readRDS)
gene_sets_list <- readRDS(here::here("6-differential_expression_analysis/tmp/gene_sets_gsea_analysis.rds"))


# Calculate signatures
signatures <- list(
  oxphos = gene_sets_list$`GO:0006119`,
  bcr_signaling = gene_sets_list$`GO:0050853`
)
mt_genes <- str_c("MT-", signatures$oxphos, sep = "") %in% rownames(seurat_list$`19`)
signatures$oxphos[mt_genes] <- str_c("MT-", signatures$oxphos[mt_genes], sep = "")
seurat_list <- purrr::map(seurat_list, function(seurat_obj) {
  seurat_obj <- AddModuleScore(seurat_obj, features = signatures, name = names(signatures))
  seurat_obj
})


# Project signatures UMAP
umaps <- purrr::map(seurat_list, function(seurat_obj) {
  umap_oxphos <- FeaturePlot(seurat_obj, features = "oxphos1", pt.size = 0.0125)
  umap_oxphos <- umap_oxphos +
    ggtitle("OXPHOS Score") +
    scale_color_viridis_c(option = "magma") +
    theme(
      plot.title = element_text(face = "plain", size = 9),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      # legend.position = c(0, 0.75),
      legend.position = "right",
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.25, "cm")
    )
  umap_bcr <- FeaturePlot(seurat_obj, features = "bcr_signaling2", pt.size = 0.0125)
  umap_bcr <- umap_bcr +
    ggtitle("BCR Signaling Score") +
    scale_color_viridis_c(option = "magma") +
    theme(
      plot.title = element_text(face = "plain", size = 9),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      # legend.position = c(0, 0.75),
      legend.position = "right",
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.3, "cm")
    )
  umaps <- (umap_oxphos / umap_bcr)
  umaps
})


# Ridge plots
ridgeplots_oxphos <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, x) {
  ridgeplot_oxphos <- plot_ridge(
    seurat_obj,
    feature = "oxphos1",
    x_lab = "OXPHOS Score",
    colors_reference = color_annotations,
    patient_id = x
  )
  ridgeplot_oxphos <- ridgeplot_oxphos +
    theme(
      axis.text = element_text(size = 5.5),
      axis.title = element_text(size = 8)
    )
  ridgeplot_oxphos
})
ridgeplots_bcr <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, x) {
  ridgeplot_bcr <- plot_ridge(
    seurat_obj,
    feature = "bcr_signaling2",
    x_lab = "BCR Signaling Score",
    colors_reference = color_annotations,
    patient_id = x
  )
  ridgeplot_bcr <- ridgeplot_bcr +
    theme(
      axis.text = element_text(size = 5.5),
      axis.title = element_text(size = 8)
    )
  ridgeplot_bcr
})


# Violin plots
## Define RT seed cells
seurat_list <- purrr::map(seurat_list, function(seurat_obj) {
  seurat_obj$is_rt <- ifelse(str_detect(seurat_obj$annotation_final, "RT"), "RT", "CLL")
  # seurat_obj$sample_description_RM <- str_c(
  #   seurat_obj$time_point,
  #   seurat_obj$sample_description,
  #   sep = "_"
  # )
  seurat_obj
})
# seurat_list$`63`$sample_description_RM <- str_c(
#   seurat_list$`63`$time_point,
#   seurat_list$`63`$sample_description_FN,
#   sep = "_"
# )


## Plot violin
color_violin <- c("#dcdddc", "#cd8899")
violins_oxphos <- purrr::map(seurat_list, function(seurat_obj) {
  violin_oxphos <- seurat_obj@meta.data %>%
    ggplot(aes(time_point, oxphos1, fill = is_rt, color = is_rt)) +
      geom_violin(alpha = 0.35, position = "dodge") +
      ggforce::geom_sina(size = 0.5, alpha = 1) +
      scale_color_manual(values = color_violin) +
      scale_fill_manual(values = color_violin) +
      labs(x = "", y = "OXPHOS\nScore", color = "", fill = "") +
      theme_classic() +
      theme(
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 8),
        legend.position = "top"
        # legend.text = element_text(size = 7)
      )
  violin_oxphos <- violin_oxphos + theme(axis.text.x = element_blank())
  violin_oxphos
})
legend_violins <- as_ggplot(get_legend(violins_oxphos$`19`))
violins_oxphos <- purrr::map(violins_oxphos, function(p) {
  p <- p + theme(legend.position = "none")
  p
})

violins_bcr <- purrr::map(seurat_list, function(seurat_obj) {
  violin_bcr <- seurat_obj@meta.data %>%
    ggplot(aes(time_point, bcr_signaling2, fill = is_rt, color = is_rt)) +
      geom_violin(alpha = 0.35) +
      ggforce::geom_sina(size = 0.5, alpha = 1) +
      scale_color_manual(values = color_violin) +
      scale_fill_manual(values = color_violin) +
      labs(x = "", y = "BCR Signaling\nScore", color = "", fill = "") +
      theme_classic() +
      theme(
        # axis.text.x = element_text(angle = 45, hjust = 0.95, color = "black", size = 7)
        axis.text.x = element_text(color = "black", size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        legend.position = "none"
        # legend.text = element_text(size = 7)
      )
  violin_bcr
})


# Arrange figure
figs <- purrr::map(names(seurat_list), function(x) {
  p <- 
    umaps[[x]] |
    (ridgeplots_oxphos[[x]] / ridgeplots_bcr[[x]]) |
    (violins_oxphos[[x]] / violins_bcr[[x]])
  p
})
names(figs) <- names(seurat_list)
fig <-
  figs$`19` / 
  figs$`365` /
  figs$`3299` /
  figs$`63`


# Save
ggsave(
  filename = here::here("results/plots/paper/rt_oxphos_and_bcr_scRNA_supplementary.pdf"),
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 29, 
  units = "cm"
)
ggsave(
  filename = here::here("results/plots/paper/rt_oxphos_and_bcr_scRNA_supplementary_legend_violins.pdf"),
  plot = legend_violins,
  device = cairo_pdf,
  width = 4, 
  height = 1, 
  units = "cm"
)