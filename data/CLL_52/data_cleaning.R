library(tidyverse)


lims_info <- read_tsv("~/Desktop/PhD/RICHTER/data/CLL_52/lims_info.txt")



lims_info$index
lims_info$flowcell
"/scratch/project/production/fastq"
"{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)


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

head(metadata_cll_52)


symbolic_links_df <- read_csv("~/Desktop/PhD/RICHTER/data/CLL_52/symlink_output_parsed.txt", col_names = FALSE)

head(symbolic_links_df)


symbolic_links_df <- symbolic_links_df[4:nrow(symbolic_links_df), ]
colnames(symbolic_links_df) <- c("symbolic_link_r1", "fastq_path_r1")

symbolic_links_df <- symbolic_links_df[str_detect(symbolic_links_df$symbolic_link_r1, "_R1_"), ]


dim(symbolic_links_df)
dim(metadata_cll_52)

metadata_cll_52 <- left_join(metadata_cll_52, symbolic_links_df, by = "fastq_path_r1")
metadata_cll_52$symbolic_link_r2 <- str_replace(metadata_cll_52$symbolic_link_r1, "_R1_", "_R2_")

zUMI <- read_tsv("~/Desktop/PhD/RICHTER/data/CLL_52/reads_for_zUMIs.samples.txt")
colnames(zUMI) <- c("symbolic_link_r1", "symbolic_link_r1", "sample", "zUMI_barcode")
zUMI <- zUMI[, c("symbolic_link_r1", "zUMI_barcode")]
metadata_cll_52_3 <- left_join(metadata_cll_52, zUMI, by = "symbolic_link_r1")


library_id_plate <- read_csv("~/Desktop/PhD/RICHTER/data/CLL_52/CLL_52_library_id_plate.csv")
indices1 <- which(library_id_plate$plate == "P2470")
indices2 <- which(library_id_plate$plate == "P2471")
library_id_plate$plate[indices1] <- "P2471"
library_id_plate$plate[indices2] <- "P2470"
library_id_plate$plate
write_csv(library_id_plate, "~/Desktop/PhD/RICHTER/data/CLL_52/CLL_52_library_id_plate.csv")
metadata_cll_52_4 <- left_join(metadata_cll_52_3, library_id_plate, by = "library_id")



gustavo_df <- read_tsv("~/Desktop/PhD/RICHTER/data/CLL_52/CLL_52_experimental_design.tsv", col_names = FALSE)
colnames(gustavo_df) <- c("original_barcode", "well")
# metadata_cll_52_4$

metadata_cll_52_4$original_barcode <- metadata_cll_52_4$plate %>%
  str_c(metadata_cll_52_4$index, sep = "_") %>%
  str_remove("-NX-xt")

all(metadata_cll_52_4$original_barcode %in% gustavo_df$original_barcode)
all(gustavo_df$original_barcode %in% metadata_cll_52_4$original_barcode)

metadata_cll_52_5 <- left_join(metadata_cll_52_4, gustavo_df, by = "original_barcode")

write_csv(metadata_cll_52_5, "~/Desktop/PhD/RICHTER/data/CLL_52/CLL_52_metadata.csv")


for (i in colnames(metadata_cll_52_4)) {
  print(any(is.na(metadata_cll_52_4[, i])))
}

unique(metadata_cll_52_4$plate[!(metadata_cll_52_4$original_barcode %in% gustavo_df$original_barcode)])


########################

metadata_cll_52_4$original_barcode[!(metadata_cll_52_4$original_barcode %in% gustavo_df$original_barcode)]


metadata_cll_52_4$original_barcode[!(metadata_cll_52_4$original_barcode %in% gustavo_df$original_barcode)] %>%
  str_subset("P2471") %>%
  str_remove("P2471_") %>%
  sort()

gustavo_df$original_barcode[!(gustavo_df$original_barcode %in% metadata_cll_52_4$original_barcode)] %>%
  str_subset("P2471") %>%
  str_remove("P2471_") %>%
  sort()



###############

lims_info <- as.data.frame(lims_info)
rownames(lims_info) <- lims_info$sample
lims_info <- lims_info[library_id_plate$library_id, ]


library_id_plate$index <- lims_info$index


library_id_plate$test_barcode <- str_c(library_id_plate$plate, library_id_plate$index, sep = "_")
library_id_plate$test_barcode <- str_remove(library_id_plate$test_barcode, "-NX-xt")

gustavo_df2 <- read_tsv("~/Desktop/PhD/RICHTER/data/CLL_52/GRCh38/data/cells.tsv", col_names = FALSE)
colnames(gustavo_df2) <- c("original_barcode", "well")
library_id_plate$test_barcode %in% gustavo_df2$original_barcode
gustavo_df2$original_barcode[!(gustavo_df2$original_barcode %in% library_id_plate$test_barcode)]

metadata_cll_52_4[metadata_cll_52_4$library_id == "AG0681", ]


all(lims_info$sample %in% library_id_plate$library_id)
all(library_id_plate$library_id %in% lims_info$sample)
library_id_plate$library_id


test <- metadata_cll_52_4
test$plate[str_detect(test$plate == "")]

indices1 <- which(test$plate == "P2470")
indices2 <- which(test$plate == "P2471")

test$plate[indices1] <- "P2471"
test$plate[indices2] <- "P2470"
test$test_barcode <- str_c(test$plate, test$index, sep = "_")
test$test_barcode <- str_remove(test$test_barcode, "-NX-xt")

all(test$test_barcode %in% gustavo_df2$original_barcode)
all(gustavo_df2$original_barcode %in% test$test_barcode)


gustavo_df2$original_barcode[!(gustavo_df2$original_barcode %in% test$test_barcode)]




metadata$sample_id <- metadata$donor_id
metadata$sample_id <- str_remove(metadata$sample_id, "CLL_")
metadata$sample_id <- str_replace(metadata$sample_id, "-", "/")
metadata$donor_id <- "63"


metadata
write_csv(metadata, "~/Desktop/PhD/RICHTER/data/CLL_52/CLL_52_metadata.csv")





