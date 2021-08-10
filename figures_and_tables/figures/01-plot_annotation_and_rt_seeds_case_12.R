# This script plots the UMAPs of the annotation and RT seed cells for case 12

# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)


# Source scripts and variables
source(here::here("bin/utils.R"))


# Load Seurat object
path_to_12 <- here::here("results/R_objects/6.seurat_annotated_12.rds")
seurat <- readRDS(path_to_12)


# Individual panels
umap_annotation <- plot_annotation(
  seurat_obj = seurat,
  pt_size = 0.5,
  colors_reference = color_annotations,
  patient_id = "12",
  nothing = TRUE
) 
umap_annotation <- umap_annotation +
  theme(
    legend.position = c(0, 0.75),
    legend.text = element_text(size = 7)
  )
genes_interest <- c("CXCR4", "CD24", "CD27", "MIR155HG", "CCND2", "PCNA",
                    "MKI67", "MZB1", "IGHM", "XBP1")
dot_plot <- plot_dot_plot(
  seurat,
  goi = rev(genes_interest),
  colors_reference = color_annotations,
  patient_id = "12"
)
dot_plot <- dot_plot +
  scale_size_continuous(range = c(0.1, 4.5)) +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = "right"
    # legend.box = "vertical"
  )
violin_plot_s_phase <- plot_violin_plot(
  seurat_obj = seurat,
  continuous_var = "S.Score",
  ylab = "S Phase Score",
  colors_reference = color_annotations,
  patient_id = "12"
)
violin_plot_s_phase <- violin_plot_s_phase +
  theme(
    axis.text.x = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25),
    axis.title.y = element_text(size = 9)
  )
violin_plot_g2m_phase <- plot_violin_plot(
  seurat_obj = seurat,
  continuous_var = "G2M.Score",
  ylab = "G2M Phase Score",
  colors_reference = color_annotations,
  patient_id = "12"
)
violin_plot_g2m_phase <- violin_plot_g2m_phase +
  theme(
    axis.title.y = element_text(size = 9),
    axis.text.x = element_blank(),
    axis.line = element_line(size = 0.25),
    axis.ticks = element_line(size = 0.25)
  )
umaps_seed_cells <- plot_split_annotation(
  seurat,
  pt_size = 0.8,
  split_by = "time_point",
  colors_reference = color_annotations,
  patient_id = "12",
  n_col = 2
)
umaps_seed_cells <- umaps_seed_cells &
  NoLegend() &
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
  ) &
  labs(x = "UMAP1", y = "UMAP2") &
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5))

fig_right <- ((violin_plot_s_phase / violin_plot_g2m_phase) | umaps_seed_cells)
fig_right <- fig_right + plot_layout(widths = c(0.5, 1))
fig_mid_and_right <- (dot_plot | fig_right) + plot_layout(widths = c(0.5, 1))
fig <- (umap_annotation | fig_mid_and_right) + plot_layout(widths = c(0.25, 0.75))
fig



# fig_left <- (umap_annotation | dot_plot)
# fig_right <- ((violin_plot_s_phase / violin_plot_g2m_phase) | umaps_seed_cells)
# fig_right <- fig_right + plot_layout(widths = c(1, 1.5))
# fig <- fig_left | fig_right

# fig_top <- (umap_annotation | dot_plot)
# fig_top[[1]] <- fig_top[[1]] +
#   theme(axis.title.x = element_text(margin = margin(t = -25, unit = "pt")))
# fig_top <- fig_top + plot_layout(widths = c(1.75, 1))
# fig_bottom <- (violin_plot_s_phase / violin_plot_g2m_phase) | umaps_seed_cells
# fig_bottom[[2]] <- fig_bottom[[2]] +
#   theme(axis.title.x = element_text(margin = margin(t = -50, unit = "pt")))
# fig_bottom <- fig_bottom + plot_layout(widths = c(1, 2))
# fig <- fig_top / fig_bottom


# To complete it, we need the total number of cells and percentage of RT cells
# per time-point
n_cells_rt <- seurat@meta.data %>%
  group_by(time_point, annotation_final) %>%
  dplyr::count(.drop = FALSE, ) %>%
  ungroup() %>%
  group_by(time_point) %>%
  mutate(
    n_cells_total = sum(n),
    percentage_total = n / n_cells_total * 100
  )
pct_seed_cells <- n_cells_rt %>%
  mutate(is_richter = str_detect(annotation_final, "RT")) %>%
  group_by(is_richter, time_point) %>%
  summarise(n_rt_cells = sum(n), pct_rt_cells = sum(percentage_total)) %>%
  filter(is_richter)
pct_seed_cells$total_n_cells_time_point <- n_cells_rt %>%
  group_by(time_point) %>%
  top_n(1) %>%
  pull(n_cells_total)
selected_cols <- c("time_point", "total_n_cells_time_point", "n_rt_cells",
                   "pct_rt_cells")
pct_seed_cells <- pct_seed_cells[, selected_cols] 


# Save
ggsave(
  filename = here::here("results/plots/paper/01-rt_annotation_and_seeds.pdf"),
  plot = fig,
  device = cairo_pdf,
  width = 21, 
  height = 11, 
  units = "cm"
)
write_delim(
  pct_seed_cells,
  file = here::here("results/plots/paper/rt_annotation_and_seeds_number_of_cells.csv")
)
