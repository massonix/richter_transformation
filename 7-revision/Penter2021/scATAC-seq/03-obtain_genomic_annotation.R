# To obtain genomic annotations, the function seqlevelsStyle needs
# access to the Internet. Because the jobs running in the cluster
# do not have access, we will run this locally to save the annotation
# and load it in future scripts


# Load packages
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)


# Obtain genomic annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"


# Save
dir.create(here::here("7-revision/Penter2021/scATAC-seq/tmp/"))
saveRDS(
  annotation,
  here::here("7-revision/Penter2021/scATAC-seq/tmp/genomic_annotation.rds")
)