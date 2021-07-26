# This script plots the UMAPs of the annotation and RT seed cells for all
# cases except 12


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)


# Source scripts and variables
source(here::here("bin/utils.R"))


# Load Seurat object
path_to_19 <- here::here("results/R_objects/6.seurat_annotated_19.rds")
path_to_63 <- here::here("results/R_objects/patient_63/3.seurat_annotated.rds")
path_to_365 <- here::here("results/R_objects/6.seurat_annotated_365.rds")
path_to_3299 <- here::here("results/R_objects/6.seurat_annotated_3299.rds")
paths_to_load <- c(
  "19" = path_to_19,
  "63" = path_to_63,
  "365" = path_to_365,
  "3299" = path_to_3299
)
seurat_list <- purrr::map(paths_to_load, readRDS)


# UMAP tissue
color_tissue <- c(
  "PB" = "#FFDBB8",
  "LN" = "#CAF9E2",
  "BM" = "#bfe10e"
)
umaps_tissue <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, x) {
  Idents(seurat_obj) <- "tissue"
  p <- DimPlot(seurat_obj, pt.size = 0.5)
  p <- p +
    scale_color_manual(
      values = color_tissue[rev(unique(seurat_obj$tissue))],
      breaks = rev(unique(seurat_obj$tissue))
    )
  p <- p +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "top",
      legend.justification = "center",
      legend.text = element_text(size = 7)
    )
  p
})
umaps_tissue <- purrr::map(umaps_tissue, function(p) {
  p <- p + theme(legend.position = "none")
  p
})
                            

# UMAPs and dot plot annotations
umaps_annotation <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, x) {
  p <- plot_annotation(
    seurat_obj = seurat_obj,
    pt_size = 0.5,
    colors_reference = color_annotations,
    patient_id = x,
    nothing = TRUE
  ) 
  p <- p +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 7)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  p
})
umaps_annotation <- purrr::map(umaps_annotation, function(p) {
  p <- p + theme(legend.position = "none")
  p
})
dot_plots <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, x) {
  p <- plot_dot_plot(
    seurat_obj,
    goi = rev(genes_to_dotplot[[x]]),
    colors_reference = color_annotations,
    patient_id = x
  )
  p <- p +
    scale_size_continuous(range = c(0.1, 4.5)) +
    theme(
      axis.text.x = element_text(size = 6),
      axis.text.y = element_blank(),
      axis.line = element_line(size = 0.25),
      axis.ticks = element_line(size = 0.25),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.position = "right",
      legend.box = 'vertical'
      # legend.box = "vertical"
    )
  p
})
dot_plots <- purrr::map(dot_plots, function(p) {
  p <- p + theme(legend.position = "none")
  p
})
dot_plots$`19` <- dot_plots$`19` +
  scale_size_continuous(range = c(0.1, 3.5))
dot_plots$`3299` <- dot_plots$`3299` +
  scale_size_continuous(range = c(0.1, 3.5))


# Violin plots proliferation
violin_plots_s_phase <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, x) {
  p <- plot_violin_plot(
    seurat_obj = seurat_obj,
    continuous_var = "S.Score",
    ylab = "S Phase\nScore",
    colors_reference = color_annotations,
    patient_id = x
  )
  p <- p +
    theme(
      axis.text.x = element_blank(),
      axis.line = element_line(size = 0.25),
      axis.ticks = element_line(size = 0.25),
      axis.title.y = element_text(size = 8)
    )
  p
})
violin_plots_g2m_phase <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, x) {
  p <- plot_violin_plot(
    seurat_obj = seurat_obj,
    continuous_var = "G2M.Score",
    ylab = "G2M Phase\nScore",
    colors_reference = color_annotations,
    patient_id = x
  )
  p <- p +
    theme(
      axis.title.y = element_text(size = 8),
      axis.text.x = element_blank(),
      axis.line = element_line(size = 0.25),
      axis.ticks = element_line(size = 0.25)
    )
  p
})


# Arrange figure
fig_rows <- purrr::map(names(seurat_list), function(x) {
  violins <- (violin_plots_s_phase[[x]] / violin_plots_g2m_phase[[x]])
  ps <- umaps_annotation[[x]] + umaps_tissue[[x]] + dot_plots[[x]] + violins
  ps <- ps +
    plot_layout(ncol = 4, widths = c(2, 1.5, 1.25, 1.25))
  ps
})
names(fig_rows) <- names(seurat_list)
fig <- (
  fig_rows$`19` /
  fig_rows$`63` /
  fig_rows$`365` /
  fig_rows$`3299` 
)


# Save
ggsave(
  filename = here::here("results/plots/paper/rt_annotation_supplementary.pdf"),
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 28, 
  units = "cm"
)


