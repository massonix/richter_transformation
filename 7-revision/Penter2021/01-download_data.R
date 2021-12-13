# Load packages
library(GEOquery)
library(glue)


# Download data
directory <- "~/Desktop/richter_transformation/data"
gse <- getGEO("GSE165087", GSEMatrix = FALSE)
gsm_list <- GSMList(gse)
for (gsm_name in names(gsm_list)) {
  print(gsm_name)
  download_dir <- glue("{directory}/{gsm_name}")
  dir.create(download_dir)
  gsm <- gsm_list[[gsm_name]]
  for (i in 1:3) {
    url <- gsm@header[[glue("supplementary_file_{i}")]]
    command <- glue::glue("wget -P {download_dir} {url}")
    system(
      command,
      intern = FALSE,
      ignore.stdout = TRUE,
      ignore.stderr = TRUE,
      wait = TRUE
    )
  }
}

