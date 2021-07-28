# This script plots the results of the differential expression analysis (RT vs CLL)


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)
library(ggrepel)
library(gridGraphics)


# Source scripts and variables
source(here::here("bin/utils.R"))


# Load data
path_to_12 <- here::here("results/R_objects/6.seurat_annotated_12.rds")
path_to_19 <- here::here("results/R_objects/6.seurat_annotated_19.rds")
path_to_63 <- here::here("results/R_objects/patient_63/3.seurat_annotated.rds")
path_to_365 <- here::here("results/R_objects/6.seurat_annotated_365.rds")
path_to_3299 <- here::here("results/R_objects/6.seurat_annotated_3299.rds")
paths_to_load <- c(
  "12" = path_to_12,
  "19" = path_to_19,
  "3299" = path_to_3299,
  "365" = path_to_365,
  "63" = path_to_63
)
seurat_list <- purrr::map(paths_to_load, readRDS)
path_to_dea_results <- here::here("results/tables/DEA/patient_specific_differential_expression_analysis_rt_vs_cll.xlsx")
dea_list <- purrr::map(names(paths_to_load), function(x) {
  df <- openxlsx::read.xlsx(xlsxFile = path_to_dea_results, sheet = x)
})
names(dea_list) <- names(paths_to_load)


# MA plots
alpha <- 0.05
min_avg_log2FC <- 0.25
selected_avg_log2FC_l <- c(
  "12" = 1.25,
  "19" = 0.75,
  "63" = 1.25,
  "365" = 1.5,
  "3299" = 0.75
)
selected_avg_log2FC <- 1.25
selected_pct_cells <- 10
selected_significance_alpha <- alpha
text_size <- 1.25
ma_plots <- purrr::map2(dea_list, names(dea_list), function(df, x) {
  selected_avg_log2FC <- selected_avg_log2FC_l[x]
  p <- ma_plot(
    df,
    selected_avg_log2FC = selected_avg_log2FC,
    selected_pct_cells = selected_pct_cells,
    selected_significance_alpha = selected_significance_alpha,
    text_size = text_size
  ) +
    ggtitle(x) +
    ylab(bquote("Mean"~log[2]~"(RT / CLL)")) +
    theme(plot.title = element_text(hjust = 0.5))
  p
})
names(ma_plots) <- names(dea_list)


# Upset plots
upregulated_richter <- purrr::map(dea_list, function(df) df$gene[df$direction == "up"])
downregulated_richter <- purrr::map(dea_list, function(df) df$gene[df$direction == "down"])
upset_upregulated <- upset(fromList(upregulated_richter), order.by = "freq")
upset_downregulated <- upset(fromList(downregulated_richter), order.by = "freq")


# Arrange figure
fig_ma_plots <-
  ma_plots$`12` + ma_plots$`19` +
  ma_plots$`63` + ma_plots$`365` +
  ma_plots$`3299` +
  plot_layout(ncol = 3)
fig_ma_plots <- fig_ma_plots &
  theme(
    legend.position = "none",
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    plot.title = element_text(size = 9)
  )


# fig <- ma_plots$`12` + ma_plots$`19` + ma_plots$`63` + ma_plots$`365` + ma_plots$`3299`
# fig <- fig + ~upset_upregulated + ~upset_downregulated
# fig <- fig + plot_layout(ncol = 3, byrow = TRUE)


# Save
ggsave(
  filename = here::here("results/plots/paper/rt_dea_ma_plots.pdf"),
  plot = fig_ma_plots,
  device = cairo_pdf,
  width = 21, 
  height = 14, 
  units = "cm"
)



