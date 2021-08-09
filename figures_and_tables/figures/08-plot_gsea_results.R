# This script plots the results of the gene set enrichment analysis (GSEA)


# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggtext)
library(ggrepel)
library(ggpubr)
library(gridGraphics)
library(clusterProfiler)
library(UpSetR)


# Source scripts and variables
source(here::here("bin/utils.R"))
# source("~/Desktop/CNAG/mnt_clust/RICHTER/current/bin/utils.R")


# Load data
path_to_gsea <- here::here("6-differential_expression_analysis/tmp/patient_specific_gsea_rt_vs_cll.rds")
path_to_gsea_raw <- here::here("6-differential_expression_analysis/tmp/patient_specific_gsea_raw_rt_vs_cll.rds")
# path_to_gsea <- "~/Desktop/CNAG/mnt_clust/RICHTER/current/6-differential_expression_analysis/tmp/patient_specific_gsea_rt_vs_cll.rds"
# path_to_gsea_raw <- "~/Desktop/CNAG/mnt_clust/RICHTER/current/6-differential_expression_analysis/tmp/patient_specific_gsea_raw_rt_vs_cll.rds"
gsea_sorted <- readRDS(path_to_gsea)
gsea_patient_specific <- readRDS(path_to_gsea_raw)


# Upset plots
upregulated_terms <- purrr::map(gsea_sorted, function(df) {
  x <- df$ID[df$NES > 0]
  x
})
downregulated_terms <- purrr::map(gsea_sorted, function(df) {
  x <- df$ID[df$NES < 0]
  x
})
upset_upregulated <- upset(
  fromList(upregulated_terms),
  order.by = "freq",
  mb.ratio = c(0.5, 0.5),
  show.numbers = FALSE
)
upset_downregulated <- upset(
  fromList(downregulated_terms),
  mb.ratio = c(0.5, 0.5),
  order.by = "freq",
  show.numbers = FALSE
)


# GSEA plots
selected_terms <- c("GO:0006119", "GO:0032543", "GO:0050853")
selected_titles <- c("oxidative phosphorylation", "mitochondrial translation",
                     "B cell receptor signaling pathway")
gsea_patient_specific <- gsea_patient_specific[c("12", "63", "365")]
gsea_plots <- purrr::map(gsea_patient_specific, function(obj) {
  plots <- purrr::map2(selected_terms, selected_titles, function(x, title) {
    p <- gseaplot(obj, geneSetID = x, by = "runningScore", title = title)
    p <- p +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.margin = margin(l = 0.25, t = 0.3, b = 0.5, r = 0.25, unit = "cm")
      )
  })
  plots
})


# Arrange figure
row1_oxphos <- gsea_plots$`12`[[1]] | gsea_plots$`63`[[1]] | gsea_plots$`365`[[1]]
row2_mt_translation <- gsea_plots$`12`[[2]] | gsea_plots$`63`[[2]] | gsea_plots$`365`[[2]] 
row3_bcr <- gsea_plots$`12`[[3]] | gsea_plots$`63`[[3]] | gsea_plots$`365`[[3]] 
fig_gsea <- row1_oxphos / row2_mt_translation / row3_bcr


# Save
ggsave(
  filename = here::here("results/plots/paper/08.1-rt_gsea_enrichment_plots.pdf"),
  plot = fig_gsea,
  device = cairo_pdf,
  width = 20, 
  height = 13.5,
  units = "cm"
)
pdf(
  file = here::here("results/plots/paper/08.2-rt_gsea_upregulated_upset_plots_supplementary.pdf"),
  onefile = FALSE,
  width = 3.94,
  height = 3
)
upset_upregulated
dev.off()
pdf(
  file = here::here("results/plots/paper/08.3-rt_gsea_downregulated_upset_plots_supplementary.pdf"),
  onefile = FALSE,
  width = 3.94,
  height = 3
)
upset_downregulated
dev.off()


# NES and q
gsea_sorted <- gsea_sorted[c("12", "63", "365")]
nes_q_dfs <- purrr::map(gsea_sorted, function(df) {
  df_sub <- df %>%
    dplyr::filter(ID %in% selected_terms) %>%
    dplyr::select("ID", "Description", "NES", "qvalues")
  df_sub
})
nes_q_df <- bind_rows(nes_q_dfs, .id = "case")
write_delim(
  nes_q_df,
  file = here::here("results/plots/paper/08-rt_gsea_nes_and_qvalues.csv"),
  delim = ";"
)
