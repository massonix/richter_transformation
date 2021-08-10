# This script plots the supplementary figure with the results of the time point-specific
# clustering


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
# source("~/Desktop/CNAG/mnt_clust/RICHTER/current/bin/utils.R")


# More functions
plot_annotation2 <- function(seurat_obj,
                             pt_size = 0.5,
                             colors_reference,
                             patient_id,
                             cells_highlight,
                             pt_size_highlight,
                             col_highlight,
                             nothing = TRUE) {
  Idents(seurat_obj) <- "annotation_final"
  p <- DimPlot(
    seurat_obj,
    pt.size = pt_size,
    order = rev(colors_reference[[patient_id]][["annotation_final"]]),
    cells.highlight = cells_highlight,
    sizes.highlight = pt_size_highlight,
    cols.highlight = col_highlight
  )
  
  if (nothing == TRUE) {
    p <- p +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
    p
  } else{
    p
  }
}


# Load data
path_to_obj <- here::here("results/R_objects/7.seurat_time_point_specific_list.rds")
# path_to_obj <- "~/Desktop/CNAG/mnt_clust/RICHTER/current/results/R_objects/7.seurat_time_point_specific_list.rds"
seurat_list_list <- readRDS(path_to_obj)


# Clean
names(seurat_list_list$`12`) <- str_remove(names(seurat_list_list$`12`), "_PB")
seurat_list_list$`12`$T6 <- NULL
names(seurat_list_list$`19`) <- str_remove(names(seurat_list_list$`19`), "_PB")
seurat_list_list$`19`$T6 <- NULL
seurat_list_list$`63`$T3_LN <- NULL
seurat_list_list$`365`$T3_LN <- NULL
names(seurat_list_list$`365`) <- str_remove(names(seurat_list_list$`365`), "_PB")
names(seurat_list_list$`3299`) <- str_remove(names(seurat_list_list$`3299`), "_BM")
seurat_list_list$`3299`$T3 <- NULL
seurat_list_list$`3299`$T2 <- NULL


# Plot
umaps_time_points <- purrr::map2(seurat_list_list, names(seurat_list_list), function(seurat_obj_list, patient_id) {
  purrr::map2(seurat_obj_list, names(seurat_obj_list), function(seurat_obj, x) {
    seurat_obj$is_rt <- ifelse(
      str_detect(seurat_obj$annotation_final, "RT"),
      "RT",
      "CLL"
    )
    p <- plot_annotation2(
      seurat_obj = seurat_obj,
      pt_size = 0.5,
      colors_reference = color_annotations,
      patient_id = patient_id,
      cells_highlight = colnames(seurat_obj)[seurat_obj$is_rt == "RT"],
      pt_size_highlight = 2,
      col_highlight = "#cd8599",
      nothing = TRUE
    )
    p$layers[[1]]$aes_params$alpha <- 0.65
    p +
      ggtitle(x) +
      theme(
        plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
        legend.position = "none"
      )
  })
})


# Arrange figure
empty_plot <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_blank() +
  ggmap::theme_nothing()
fig_row1 <- 
  (umaps_time_points$`12`$T1 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T1 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`63`$T1_PB + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (umaps_time_points$`365`$T2 + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  umaps_time_points$`3299`$T1

fig_row2 <- 
  (umaps_time_points$`12`$T2 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T3 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`63`$T1_LN + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  empty_plot


fig_row3 <- 
  (umaps_time_points$`12`$T4 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T4 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`63`$T2_PB + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  empty_plot

fig_row4 <- 
  (umaps_time_points$`12`$T5 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T5 + theme(plot.margin = margin(l = 0.25, t = 0.3, r = 1.5, unit = "cm"))) |
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  empty_plot


fig <- (
  plot_spacer() /
    fig_row1 /
    fig_row2 /
    fig_row3 /
    fig_row4 /
    plot_spacer()
)
fig <- fig + plot_layout(heights = c(0.25, 1, 1, 1, 1, 0.25))


# Save
ggsave(
  filename = here::here("results/plots/paper/07-rt_time_point_specific_clustering_supplementary.pdf"),
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 19, 
  units = "cm"
)
# ggsave(
#   filename = "~/Desktop/CNAG/mnt_clust/RICHTER/current/results/plots/paper/07-rt_time_point_specific_clustering_supplementary.pdf",
#   plot = fig,
#   device = cairo_pdf,
#   width = 21, 
#   height = 23, 
#   units = "cm"
# )


# Create and save legend
seurat_obj <- seurat_list_list$`12`$T4
seurat_obj$is_rt <- ifelse(
  str_detect(seurat_obj$annotation_final, "RT"),
  "RT",
  "CLL"
)
p <- plot_annotation2(
  seurat_obj = seurat_obj,
  pt_size = 0.5,
  colors_reference = color_annotations,
  patient_id = patient_id,
  cells_highlight = colnames(seurat_obj)[seurat_obj$is_rt == "RT"],
  pt_size_highlight = 2,
  col_highlight = "#cd8599",
  nothing = TRUE
)
p$layers[[1]]$aes_params$alpha <- 0.65
leg <- as_ggplot(get_legend(p))
ggsave(
  filename = here::here("results/plots/paper/legends/07-rt_time_point_specific_clustering_supplementary_legend.pdf"),
  plot = leg,
  device = cairo_pdf,
  width = 12, 
  height = 4, 
  units = "cm"
)
# ggsave(
#   filename = "~/Desktop/CNAG/mnt_clust/RICHTER/current/results/plots/paper/legends/07-rt_time_point_specific_clustering_supplementary_legend.pdf",
#   plot = leg,
#   device = cairo_pdf,
#   width = 12, 
#   height = 4, 
#   units = "cm"
# )
