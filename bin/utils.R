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
      size = text_size,
      max.overlaps = max_overlaps
    )
}



plot_annotation <- function(seurat_obj,
                            pt_size = 0.5,
                            colors_reference,
                            patient_id,
                            nothing = TRUE) {
  Idents(seurat_obj) <- "annotation_final"
  p <- DimPlot(seurat_obj, pt.size = pt_size)
  p <- p +
    scale_color_manual(
      values = colors_reference[[patient_id]][["colors"]],
      breaks = colors_reference[[patient_id]][["annotation_final"]],
      labels = colors_reference[[patient_id]][["custom_labels"]]
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


plot_split_annotation <- function(seurat_obj,
                                  pt_size = 0.5,
                                  split_by,
                                  colors_reference,
                                  patient_id,
                                  n_col = 3) {
  p <- DimPlot(
    seurat_obj,
    pt.size = pt_size,
    split.by = split_by,
    ncol = n_col
  ) &
    scale_color_manual(
      values = colors_reference[[patient_id]][["colors"]],
      breaks = colors_reference[[patient_id]][["annotation_final"]],
      labels = colors_reference[[patient_id]][["custom_labels"]]
    )
  p
}


plot_dot_plot <- function(seurat_obj,
                          goi,
                          colors_reference,
                          patient_id) {
  p <- DotPlot(seurat_obj, features = rev(goi)) +
    scale_color_distiller(palette = "Blues", direction = 1) +
    scale_y_discrete(
      limits = rev(colors_reference[[patient_id]][["annotation_final"]]),
      labels = rev(colors_reference[[patient_id]][["custom_labels"]])
    ) +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      legend.title = element_text(size = 12)
    )
  p$guides$colour$title <- "Average\nExpression"
  p$guides$size$title <- "Percent\nExpressed"
  p
}


plot_violin_plot <- function(seurat_obj,
                             continuous_var,
                             ylab,
                             colors_reference,
                             patient_id) {
  p <- seurat_obj@meta.data %>%
    ggplot(aes_string("annotation_final", continuous_var, fill = "annotation_final")) +
    geom_hline(yintercept = 0, color = "gray35", linetype = "dotted") +
    geom_violin(color = NA) +
    labs(x = "", y = ylab) +
    scale_x_discrete(
      breaks = colors_reference[[patient_id]][["annotation_final"]],
      labels = colors_reference[[patient_id]][["custom_labels"]]) +
    scale_fill_manual(
      values = colors_reference[[patient_id]][["colors"]],
      breaks = colors_reference[[patient_id]][["annotation_final"]],
      labels = colors_reference[[patient_id]][["custom_labels"]]
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        color = "black",
        angle = 45,
        vjust = 1,
        hjust = 1,
        size = 11
      )
    )
  p
}


plot_ridge <- function(seurat_obj,
                       feature,
                       x_lab,
                       colors_reference,
                       patient_id) {
  p <- seurat_obj@meta.data %>%
    mutate(annotation_final = factor(
      annotation_final,
      levels = rev(colors_reference[[patient_id]][["annotation_final"]])
    )) %>%
    ggplot(aes_string(feature, "annotation_final", fill = "annotation_final")) +
    ggridges::geom_density_ridges() +
    scale_y_discrete(
      breaks = rev(colors_reference[[patient_id]][["annotation_final"]]),
      labels = rev(colors_reference[[patient_id]][["custom_labels"]])
    ) +
    scale_fill_manual(
      values = rev(colors_reference[[patient_id]][["colors"]]),
      breaks = rev(colors_reference[[patient_id]][["annotation_final"]]),
      labels = rev(colors_reference[[patient_id]][["custom_labels"]])
    ) +
    labs(title = "", x = x_lab, y = "") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )
  p
}


ggpreview <- function (x, w = 5, h = 5, dpi = 150, units = "in") {
  if (!units %in% c("in", "px", "cm", "mm")) {
    stop("units has to be: in, px, cm, or mm")
  }
  tmp <- tempfile(fileext = ".png")
  grDevices::png(filename = tmp, width = w, height = h, res = dpi, 
                 units = units)
  print(x)
  grDevices::dev.off()
  rstudioapi::viewer(tmp)
}

color_annotations <- list(
  "12" = list(
    annotation_final = c("CXCR4hiCD27lo", "CXCR4loCD27hi", "MIR155HGhi", "CCND2lo RT",
                         "CCND2hi RT", "RT proliferative", "MZB1hiIGHMhiXBP1hi"),
    colors = c("#cfc9c9", "#808080", "#333333", "#fcd2cf", "#d65c7a", "#813151",
               "mediumorchid3"),
    custom_labels = list(bquote(CXCR4^hi~CD27^lo), bquote(CXCR4^lo~CD27^hi), bquote(MIR155HG^hi),
                         bquote(CCND2^lo~"RT"), bquote(CCND2^hi~"RT"), bquote("RT proliferative"),
                         bquote(MZB1^hi~IGHM^hi~XBP1^hi))
  ),
  "19" = list(
    annotation_final = c("CXCR4hiCD27lo", "CXCR4loCD27hi", "MIR155HGhi CLL",
                         "MIR155HGhi RT", "TP53INP1hi RT", "RT proliferative"),
    colors = c("#cfc9c9", "#808080", "#333333", "#fcd2cf", "#d65c7a", "#813151"),
    custom_labels = list(bquote(CXCR4^hi~CD27^lo), bquote(CXCR4^lo~CD27^hi), bquote(MIR155HG^hi~"CLL"),
                         bquote(MIR155HG^hi~"RT"), bquote(TP53INP1^hi~"RT"), 
                         bquote("RT proliferative"))
  ),
  "63" = list(
    annotation_final = c("CLL", "RT quiescent", "RT proliferative"),
    colors = c("#808080", "#d65c7a", "#813151"),
    custom_labels = list(bquote("CLL"), bquote("RT quiescent"), bquote("RT proliferative"))
  ),
  "365" = list(
    annotation_final = c("CLL", "CXCR4loCD27hi RT", "MIR155HGhi RT", "RT quiescent",
                         "RT proliferative"),
    colors = c("#cfc9c9", "#808080", "#333333", "#d65c7a", "#813151"),
    custom_labels = list(bquote("CLL"), bquote(CXCR4^lo~CD27^hi~"RT"), bquote(MIR155HG^hi~"RT"),
                         bquote("RT quiescent"), bquote("RT proliferative"))
  ),
  "3299" = list(
    annotation_final = c("CXCR4hiCD27lo", "CXCR4loCD27hi", "CD83loMIR155HGhi",
                         "CD83hiMIR155HGhi", "RT"),
    colors = c("#d4d4d4", "#8c8c8c", "#5e5e5e", "#000000", "#d65c7a"),
    custom_labels = list(bquote(CXCR4^hi~CD27^lo), bquote(CXCR4^lo~CD27^hi), bquote(CD83^lo~MIR155HG^hi),
                         bquote(CD83^hi~MIR155HG^hi), bquote("RT"))
  )
)


genes_to_dotplot <- list(
  "19" = c("CXCR4", "CD27", "S100A4", "MIR155HG", "MS4A1", "ENO1",
           "CCR7", "TCL1A", "TP53INP1", "PIK3IP1", "MKI67", "PCNA"),
  "63" = c("CXCR4", "TCL1A", "BTK", "WNT3", "PCNA", "MKI67"),
  "365" = c("CXCR4", "CD27", "MIR155HG", "ENO1", "TCL1A",
            "PCNA", "MKI67"),
  "3299" =  c("CXCR4", "CD24", "CD27", "S100A4", "MS4A1",
              "MIR155HG",  "CD83", "ENO1", "CD40", "CCR7", "KLF6",
              "PIK3IP1")
)


rt_vs_cll_dict <- list(
  "12" = c(
    "CXCR4loCD27hi" = "normal",
    "CXCR4hiCD27lo" = "normal",
    "MIR155HGhi" = "normal",
    "MZB1hiIGHMhiXBP1hi" = "normal",
    "CCND2hi RT" = "malignant",
    "CCND2lo RT" = "malignant",
    "RT proliferative" = "malignant"
  ),
  "19" = c(
    "CXCR4hiCD27lo" = "normal",
    "CXCR4loCD27hi" = "normal",
    "MIR155HGhi CLL" = "normal",
    "MIR155HGhi RT" = "malignant",
    "TP53INP1hi RT" = "malignant",
    "RT proliferative" = "malignant"
  ),
  "63" = c(
    "CLL" = "normal",
    "RT quiescent" = "malignant",
    "RT proliferative" = "malignant"
  ),
  "365" = c(
    "CLL" = "normal",
    "CXCR4loCD27hi RT" = "malignant",
    "MIR155HGhi RT" = "malignant",
    "RT quiescent" = "malignant",
    "RT proliferative" = "malignant"
  ),
  "3299" = c(
    "CXCR4hiCD27lo" = "normal",
    "CXCR4loCD27hi" = "normal",
    "CD83loMIR155HGhi" = "normal",
    "CD83hiMIR155HGhi" = "normal",
    "RT" = "malignant"
  )
)


find_assay_specific_features <- function(seurat_obj,
                                         assay_var = "assay",
                                         n_features = 5000) {
  seurat_list <- SplitObject(seurat_obj, split.by = assay_var)
  seurat_list <- purrr::map(
    seurat_list,
    FindVariableFeatures,
    nfeatures = n_features
  )
  hvg <- purrr::map(seurat_list, VariableFeatures)
  shared_hvg <- Reduce(intersect, hvg)
  shared_hvg
}
integrate_assays <- function(seurat_obj,
                             assay_specific = TRUE,
                             assay_var = "assay",
                             shared_hvg,
                             n_dim = 30
) {
  if (assay_specific) {
    seurat_obj <- seurat_obj %>%
      ScaleData(features = shared_hvg) %>%
      RunPCA(features = shared_hvg) %>%
      RunHarmony(group.by.vars = assay_var, reduction = "pca", dims = 1:n_dim)
  } else {
    seurat_obj <- seurat_obj %>%
      ScaleData() %>%
      RunPCA()
  }
  
  seurat_obj
}

