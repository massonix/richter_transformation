# This script plots the UMAPs of early RT seeding cells for all
# cases except 12


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)


# Source scripts and variables
source(here::here("bin/utils.R"))


# Parameters
x_limits <- list(
  "19" = c(-6, 9),
  "63" = c(-3.5, 3.5),
  "365" = c(-14, 6),
  "3299" = c(-8, 6.5)
)
y_limits <- list(
  "19" = c(-5, 5),
  "63" = c(-5.5, 4.5),
  "365" = c(-5, 9),
  "3299" = c(-5.1, 5.5)
)


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


# UMAP RT seeds
umaps_seed_cells <- purrr::map2(seurat_list, names(seurat_list), function(seurat_obj, patient_id) {
  df <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  df$time_point <- seurat_obj$time_point
  df$annotation_final <- seurat_obj$annotation_final
  time_points <- sort(unique(seurat_obj$time_point))
  ps <- purrr::map(time_points, function(x) {
    p <- df %>%
      filter(time_point == x) %>%
      ggplot(aes(UMAP_1, UMAP_2, color = annotation_final)) +
        geom_point(size = 0.4) +
        scale_x_continuous(limits = x_limits[[patient_id]]) +
        scale_y_continuous(limits = y_limits[[patient_id]]) +
        scale_color_manual(
          values = color_annotations[[patient_id]][["colors"]],
          breaks = color_annotations[[patient_id]][["annotation_final"]],
          labels = color_annotations[[patient_id]][["custom_labels"]]
        ) +
        ggmap::theme_nothing() +
        ggtitle(x) +
        theme(plot.title = element_text(size = 9, hjust = 0.5))
  })
  names(ps) <- time_points
  ps
})


# For case 63 we have two tissues for the same time-point, so we will plot them
# separately
df <- as.data.frame(Embeddings(seurat_list$`63`, reduction = "umap"))
df$time_point <- seurat_list$`63`$time_point
df$tissue <- seurat_list$`63`$tissue
df$time_point2 <- str_c(df$time_point, df$tissue, sep = "_")
df$annotation_final <- seurat_list$`63`$annotation_final
time_points <- sort(unique(df$time_point2))
ps <- purrr::map(time_points, function(x) {
  p <- df %>%
    filter(time_point2 == x) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = annotation_final)) +
    geom_point(size = 0.4) +
    scale_x_continuous(limits = x_limits[["63"]]) +
    scale_y_continuous(limits = y_limits[["63"]]) +
    scale_color_manual(
      values = color_annotations[["63"]][["colors"]],
      breaks = color_annotations[["63"]][["annotation_final"]],
      labels = color_annotations[["63"]][["custom_labels"]]
    ) +
    ggmap::theme_nothing() +
    ggtitle(x) +
    theme(plot.title = element_text(size = 9, hjust = 0.5))
})
names(ps) <- time_points
umaps_seed_cells$`63` <- ps


# Arrange figure
empty_plot <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_blank() +
  ggmap::theme_nothing()
fig_row1 <- 
  (umaps_seed_cells$`19`$T1 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_seed_cells$`63`$T1_PB + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (umaps_seed_cells$`365`$T2 + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  umaps_seed_cells$`3299`$T1

fig_row2 <- 
  (umaps_seed_cells$`19`$T3 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_seed_cells$`63`$T1_LN + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (umaps_seed_cells$`365`$T3 + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  umaps_seed_cells$`3299`$T2


fig_row3 <- 
  (umaps_seed_cells$`19`$T4 + theme(plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  (umaps_seed_cells$`63`$T2_PB + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  umaps_seed_cells$`3299`$T3

fig_row4 <- 
  (umaps_seed_cells$`19`$T5 + theme(plot.margin = margin(l = 0.25, t = 0.3, r = 1.5, unit = "cm"))) |
  (umaps_seed_cells$`63`$T3_LN + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm")))|
  (empty_plot + theme(plot.margin = margin(t = 0.3, b = 0.5, r = 1.5, unit = "cm"))) |
  empty_plot

fig_row5 <- 
  (umaps_seed_cells$`19`$T6 + theme(plot.margin = margin(l = 0.25, t = 0.3, r = 1.5, b = 0.5,  unit = "cm"))) |
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
  filename = here::here("results/plots/paper/04-rt_seed_cells_supplementary.pdf"),
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 23, 
  units = "cm"
)


# To complete it, we need the total number of cells and percentage of RT cells
# per time-point
pct_seed_cells_dfs <- purrr::map(seurat_list, function(seurat_obj) {
  seurat_obj@meta.data$time_point2 <- str_c(seurat_obj$time_point, seurat_obj$tissue, sep = "_")
  n_cells_rt <- seurat_obj@meta.data %>%
    group_by(time_point2, annotation_final) %>%
    dplyr::count(.drop = FALSE, ) %>%
    ungroup() %>%
    group_by(time_point2) %>%
    mutate(
      n_cells_total = sum(n),
      percentage_total = n / n_cells_total * 100
    )
  pct_seed_cells <- n_cells_rt %>%
    mutate(is_richter = str_detect(annotation_final, "RT")) %>%
    group_by(is_richter, time_point2) %>%
    summarise(n_rt_cells = sum(n), pct_rt_cells = sum(percentage_total)) %>%
    filter(is_richter)
  pct_seed_cells$total_n_cells_time_point <- n_cells_rt %>%
    group_by(time_point2) %>%
    top_n(1) %>%
    pull(n_cells_total)
  selected_cols <- c("time_point2", "total_n_cells_time_point", "n_rt_cells",
                     "pct_rt_cells")
  pct_seed_cells <- pct_seed_cells[, selected_cols]
  pct_seed_cells
})

pct_seed_cells_df <- bind_rows(pct_seed_cells_dfs, .id = "case")
write_delim(
  pct_seed_cells_df,
  file = here::here("results/plots/paper/04-rt_annotation_and_seeds_number_of_cells_supplementary.csv"),
  delim = ";"
)
