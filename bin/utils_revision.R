# This script contains all the functions needed to plot the figures
# after the revision


# Variables
color_palette <- c("black", "gray", "red", "yellow", "violet", "green4",
                   "blue", "mediumorchid2", "coral2", "blueviolet",
                   "indianred4", "deepskyblue1", "dimgray", "deeppink1",
                   "green", "lightgray", "hotpink1")
cols_rt <- c("#dcdddc", "#dc7796")


# Functions
plot_nes <- function(gsea, alpha, top_n = 5, cols = c("#8974b9", "#ff706e")) {
  df <- gsea@result %>%
    dplyr::filter(p.adjust < alpha) %>%
    dplyr::arrange(desc(NES))
  df <- df[c(1:top_n, (nrow(df) - top_n + 1):nrow(df)), ]
  p <- df %>%
    dplyr::mutate(direction = ifelse(
      NES > 0,
      "up",
      "down"
    )) %>%
    ggplot(aes(NES, fct_reorder(Description, NES), fill = direction)) +
    geom_col() +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values = cols) +
    ylab("") +
    theme_bw() +
    theme(legend.position = "none")
}


plot_signature <- function(seurat_obj, signature, group.by, fill.by, cols, ylab) {
  df <- seurat_obj@meta.data
  p <- df %>%
    ggplot(aes_string(group.by, signature, fill = fill.by)) +
    geom_violin() +
    geom_boxplot(width = 0.25, outlier.size = 0.1) +
    scale_fill_manual(values = cols) +
    ylab(ylab) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title = element_blank(),
      legend.position = "top"
    )
  p
}


score_oxphos <- function(seurat_obj, gsea, name = "oxphos_score") {
  gene_sets <- gsea@geneSets
  x <- gene_sets$`GO:0006119`
  mt_genes <- rownames(seurat) %>%
    str_subset("^MT-") %>%
    str_remove("^MT-")
  x[x %in% mt_genes] <- str_c("MT-", x[x %in% mt_genes], sep = "")
  seurat_obj <- AddModuleScore(seurat_obj, features = list(x), name = name)
  vars <- colnames(seurat_obj@meta.data)
  colnames(seurat_obj@meta.data)[vars == str_c(name, 1, sep = "")] <- name
  seurat_obj
}


score_bcr <- function(seurat_obj, gsea, name = "bcr_score") {
  gene_sets <- gsea@geneSets
  x <- gene_sets$`GO:0050853`
  seurat_obj <- AddModuleScore(seurat_obj, features = list(x), name = name)
  vars <- colnames(seurat_obj@meta.data)
  colnames(seurat_obj@meta.data)[vars == str_c(name, 1, sep = "")] <- name
  seurat_obj
}

