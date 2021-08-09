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


# Load data
path_to_obj <- here::here("results/R_objects/7.seurat_time_point_specific_list.rds")
seurat_list_list <- readRDS(path_to_obj)


# Plot
umaps_time_points <- purrr::map2(seurat_list_list, names(seurat_list_list), function(seurat_obj_list, patient_id) {
  purrr::map2(seurat_obj_list, names(seurat_obj_list), function(seurat_obj, x) {
    p <- plot_annotation(
      seurat_obj = seurat_obj,
      pt_size = 0.65,
      colors_reference = color_annotations,
      patient_id = patient_id,
      nothing = TRUE
    )
    p +
      ggtitle(x) +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "none"
      )
  })
})


# Arrange figure
empty_plot <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_blank() +
  ggmap::theme_nothing()
fig_row1 <- 
  (umaps_time_points$`12`$T1_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T1_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`63`$T1_PB + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (umaps_time_points$`365`$T2_PB + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  umaps_time_points$`3299`$T1_BM

fig_row2 <- 
  (umaps_time_points$`12`$T2_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T3_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`63`$T1_LN + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (umaps_time_points$`365`$T3_LN + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  umaps_time_points$`3299`$T2_BM


fig_row3 <- 
  (umaps_time_points$`12`$T4_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T4_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`63`$T2_PB + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  umaps_time_points$`3299`$T3_BM

fig_row4 <- 
  (umaps_time_points$`12`$T5_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T5_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`63`$T3_LN + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  empty_plot

fig_row5 <- 
  (umaps_time_points$`12`$T6_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_time_points$`19`$T6_PB + theme(plot.margin = margin(l = 0.25, t = 0.3, r = 1.5, b = 0.5,  unit = "cm"))) |
  (empty_plot + theme(plot.margin = margin(t = 0.3, r = 1.5, b = 0.5, unit = "cm")))|
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  empty_plot


fig <- (
  plot_spacer() /
    fig_row1 /
    fig_row2 /
    fig_row3 /
    fig_row4 /
    fig_row5 /
    plot_spacer()
)
fig <- fig + plot_layout(heights = c(0.25, 1, 1, 1, 1, 1, 0.25))


# Save
ggsave(
  filename = here::here("results/plots/paper/07-rt_time_point_specific_clustering_supplementary.pdf"),
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 23, 
  units = "cm"
)