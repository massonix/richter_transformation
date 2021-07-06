# This data was inherited from a 2018 project, and we need to do some 
# data archaeology to have proper metadata


# Load packages
library(tidyverse)
library(here)


# Load data
lims_info <- read_tsv(here("data/CLL_52/lims_info.txt"))
symbolic_links_df <- read_csv(
  here("data/CLL_52/symlink_output_parsed.txt"),
  col_names = FALSE
)
zUMI <- read_tsv(here("data/CLL_52/reads_for_zUMIs.samples.txt"))
library_id_plate <- read_csv(here("data/CLL_52/CLL_52_library_id_plate.csv"))
exper_design_df <- read_tsv(
  here("data/CLL_52/CLL_52_experimental_design.tsv"),
  col_names = FALSE
)


# Generate fastq paths and initialize metadata
# In Smart-seq2, every single-cell (well) is a library. Thus, each row of the 
# LIMS output corresponds to a single-cell.
path_read_1 <- str_c(
  "/scratch/project/production/fastq/",
  lims_info$flowcell,
  "/",
  lims_info$lane,
  "/fastq/",
  lims_info$flowcell,
  "_",
  lims_info$lane,
  "_",
  lims_info$index,
  "_1.fastq.gz",
  sep = ""
)
path_read_2 <- str_c(
  "/scratch/project/production/fastq/",
  lims_info$flowcell,
  "/",
  lims_info$lane,
  "/fastq/",
  lims_info$flowcell,
  "_",
  lims_info$lane,
  "_",
  lims_info$index,
  "_2.fastq.gz",
  sep = ""
)
index <- lims_info$index
library_id <- lims_info$sample
subproject <- lims_info$subproject
donor_id <- lims_info$SampleName
metadata_cll_52 <- data.frame(
  subproject = subproject,
  library_id = library_id,
  index = index,
  donor_id = donor_id,
  fastq_path_r1 = path_read_1,
  fastq_path_r2 = path_read_2
)


# Parse symbolic links
# Symlinks are needed to match the output of zUMIs to the output of LIMS,
# so that we can match the cell barcodes to the appropriate metadata
symbolic_links_df <- symbolic_links_df[4:nrow(symbolic_links_df), ]
colnames(symbolic_links_df) <- c("symbolic_link_r1", "fastq_path_r1")
symbolic_links_df <- symbolic_links_df[str_detect(symbolic_links_df$symbolic_link_r1, "_R1_"), ]
metadata_cll_52 <- left_join(metadata_cll_52, symbolic_links_df, by = "fastq_path_r1")
metadata_cll_52$symbolic_link_r2 <- str_replace(metadata_cll_52$symbolic_link_r1, "_R1_", "_R2_")


# Parse zUMIs barcodes
colnames(zUMI) <- c("symbolic_link_r1", "symbolic_link_r1", "sample", "zUMI_barcode")
zUMI <- zUMI[, c("symbolic_link_r1", "zUMI_barcode")]
metadata_cll_52 <- left_join(metadata_cll_52, zUMI, by = "symbolic_link_r1")


# Add plate info (note we interchanged plates P2470 and P2471 to match Gustavo's metadata)
metadata_cll_52 <- left_join(
  metadata_cll_52,
  library_id_plate,
  by = "library_id"
)


# Add info regarding the experimental design (well...)
colnames(exper_design_df) <- c("original_barcode", "well")
metadata_cll_52$original_barcode <- metadata_cll_52$plate %>%
  str_c(metadata_cll_52$index, sep = "_") %>%
  str_remove("-NX-xt")
all(metadata_cll_52$original_barcode %in% exper_design_df$original_barcode)
all(exper_design_df$original_barcode %in% metadata_cll_52$original_barcode)
metadata_cll_52 <- left_join(metadata_cll_52, exper_design_df, by = "original_barcode")


# Add sample_id
metadata_cll_52$sample_id <- metadata_cll_52$donor_id
metadata_cll_52$sample_id <- str_remove(metadata_cll_52$sample_id, "CLL_")
metadata_cll_52$sample_id <- str_replace(metadata_cll_52$sample_id, "-", "/")
metadata_cll_52$donor_id <- "63"


# Save
write_csv(metadata_cll_52, here("data/CLL_52/CLL_52_experiment_metadata.csv"))


