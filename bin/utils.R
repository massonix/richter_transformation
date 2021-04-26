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