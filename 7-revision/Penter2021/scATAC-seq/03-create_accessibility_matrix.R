library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
frag_path <- "~/Desktop/richter_transformation/data/scATAC/GSM4982189/GSM4982189_CLL9_1.fragments.tsv.gz"
frag_obj <- CreateFragmentObject(frag_path)
macs_path <- "/Users/rmassoni/opt/miniconda3/bin/macs2"
peaks <- CallPeaks(frag_obj, macs2.path = macs_path)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(
  x = peaks,
  ranges = blacklist_hg38_unified,
  invert = TRUE
)
mat <- FeatureMatrix(fragments = frag_obj, features = peaks)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
chrom_assay <- CreateChromatinAssay(
  counts = mat,
  ranges = peaks,
  genome = "hg38",
  fragments = frag_obj,
  annotation = annotation,
  sep = c("-", "-"),
  min.features = 200
)
seurat_obj <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")


# zcat < GSM4982189_CLL9_1.fragments.tsv.gz | cut -f4 | sort | uniq -c > fragments_per_cell.txt
# zcat < GSM4982189_CLL9_1.fragments.tsv.gz | cut -f4 | sort | uniq -c | awk '{print $2, $1}' > fragments_per_cell.txt
path_to_n_fragments <- "/Users/rmassoni/Desktop/richter_transformation/data/scATAC/GSM4982189/fragments_per_cell.txt"
n_frag_df <- read_delim(
  path_to_n_fragments,
  col_names = c("cell_barcode", "n_fragments")
)
n_frag_df <- n_frag_df %>%
  arrange(desc(n_fragments)) %>%
  mutate(position = 1:nrow(n_frag_df))
n_frag_df %>%
  ggplot(aes(position, n_fragments)) +
    geom_line() +
    labs(x = "barcodes", y = "Number of Fragments") +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    geom_vline(xintercept = 5350, linetype = "dashed", color = "red")

