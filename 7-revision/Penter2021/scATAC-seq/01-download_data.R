# This script downloads the scATAC-seq of case CLL9 from 10.1158/2159-8290
# available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163579

library(tidyverse)
library(glue)
library(GEOquery)


# Download data
# directory <- "~/Desktop/richter_transformation/data/scATAC"
directory1 <- here::here("data/Penter2021/scATAC-seq")
dir.create(directory1)
gse <- getGEO("GSE163579", GSEMatrix = FALSE)
saveRDS(gse, here::here("results/R_objects/gse_atac_obj.rds"))
gsm_list <- GSMList(gse)
supplementary_files <- purrr::map_chr(gsm_list, function(x) {
  x@header$supplementary_file_1
})
is_richter <- str_detect(supplementary_files, "CLL9")
supplementary_files <- supplementary_files[is_richter]
for (gsm_name in names(supplementary_files)) {
  print(gsm_name)
  download_dir <- glue("{directory1}/{gsm_name}")
  dir.create(download_dir, recursive = TRUE)
  url <- supplementary_files[gsm_name]
  command <- glue("wget -P {download_dir} {url}")
  system(
    command,
    intern = FALSE,
    ignore.stdout = TRUE,
    ignore.stderr = TRUE,
    wait = TRUE
  )
}