# Useful functions used in different scripts


horizontal_barplot <- function(df, categorical_var, continuous_var, ylab) {
  levels_var <- df[[categorical_var]][order(df[[continuous_var]])]
  df[[categorical_var]] <- factor(df[[categorical_var]], levels = levels_var)
  df %>%
    ggplot(aes_string(categorical_var, continuous_var)) +
      geom_col() +
      labs(x = "", y = ylab) +
      theme_bw() +
      theme(axis.text = element_text(size = 11),
            axis.title = element_text(size = 13)) +
      coord_flip()
}



customized_boxplot <- function(p) {
  p +
    geom_boxplot(width = 0.75) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 13),
          axis.title.y = element_text(size = 13))
}



horizontal_boxplot <- function(df,
                               categorical_var,
                               continuous_var,
                               fill,
                               ylab,
                               decreasing = FALSE) {
  unordered_lev <- unique(df[[categorical_var]])
  means_cont <- purrr::map_dbl(unordered_lev, function(x) {
    mean(df[[continuous_var]][df[[categorical_var]] == x])
  })
  names(means_cont) <- unordered_lev
  ordered_lev <- names(means_cont)[order(means_cont, decreasing = decreasing)]
  df[[categorical_var]] <- factor(df[[categorical_var]], levels = ordered_lev)
  df %>%
    ggplot(aes_string(categorical_var, continuous_var, fill = fill)) +
    geom_boxplot() +
    labs(x = "", y = ylab) +
    theme_bw() +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 13)) +
    coord_flip()
}



horizontal_violin_plot <- function(df,
                               categorical_var,
                               continuous_var,
                               fill,
                               ylab,
                               decreasing = FALSE) {
  unordered_lev <- unique(df[[categorical_var]])
  means_cont <- purrr::map_dbl(unordered_lev, function(x) {
    mean(df[[continuous_var]][df[[categorical_var]] == x])
  })
  names(means_cont) <- unordered_lev
  ordered_lev <- names(means_cont)[order(means_cont, decreasing = decreasing)]
  df[[categorical_var]] <- factor(df[[categorical_var]], levels = ordered_lev)
  df %>%
    ggplot(aes_string(categorical_var, continuous_var, fill = fill)) +
    geom_violin() +
    labs(x = "", y = ylab) +
    theme_bw() +
    theme(axis.text = element_text(size = 11),
          axis.title = element_text(size = 13)) +
    coord_flip()
}



plot_histogram_qc <- function(df, x, x_lab) {
  df %>%
    ggplot(aes_string(x)) +
    geom_histogram(bins = 100) +
    labs(x = x_lab, y = "Number of Cells") +
    theme_pubr()
}


process_seurat <- function(seurat_obj, dims = 1:30) {
  seurat_obj <- seurat_obj %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(reduction = "pca", dim = dims)
  seurat_obj
}


plot_split_umap <- function(seurat_obj, var, pt_size) {
  df <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
  df$new_var <- seurat_obj@meta.data[[var]]
  p <- df %>%
    ggplot(aes(UMAP_1, UMAP_2, color = new_var)) +
    geom_point(size = pt_size) +
    facet_wrap(~new_var) +
    theme_classic() +
    theme(legend.position = "none")
  p
}


run_enrichr <- function(x,
                        database = "GO_Biological_Process_2018",
                        max_total_genes,
                        min_enriched_genes,
                        p_adj_threshold,
                        min_odds_ratio) {
  go <- enrichr(genes = x, databases = database)
  go <- go[[database]]
  go <- separate(
    go,
    col = "Overlap",
    sep = "/",
    into = c("n_genes_enriched", "n_genes_total")
  )
  selected_terms <- 
    go$n_genes_total < max_total_genes &
    go$n_genes_enriched >= min_enriched_genes &
    go$Adjusted.P.value < p_adj_threshold &
    go$Odds.Ratio > min_odds_ratio
  go <- go[selected_terms, ]
  go <- arrange(go, desc(Odds.Ratio))
  go
}


ma_plot <- function(df,
                    point_size = 0.75,
                    point_alpha = 0.8,
                    cols = c("#cc3536", "#2d63a1", "gray"),
                    selected_avg_log2FC,
                    selected_pct_cells,
                    selected_significance_alpha,
                    text_size = 2.5,
                    max_overlaps = 20) {
  p <- df %>%
    mutate(direction = factor(
      direction,
      levels = c("up", "down", "not sig.")
    )) %>%
    ggplot(aes(pct_cells_expressing, avg_log2FC, color = direction)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = c("#cc3536", "#2d63a1", "gray")) +
    labs(
      x = "Percent Expressed (%)",
      y = "Mean log2(RT / CLL)",
      color = ""
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 11)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2)))
  selected_data <- df %>%
    filter(
      abs(avg_log2FC) >= selected_avg_log2FC &
        pct_cells_expressing > selected_pct_cells &
        p_val_adj < selected_significance_alpha
    )
  p +
    geom_text_repel(
      data = selected_data,
      aes(label = gene),
      color = "black",
      size = 2.5,
      max.overlaps = 20
    )
}
