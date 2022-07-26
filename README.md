# Detection of early seeding of Richter transformation in chronic lymphocytic leukemia

Richter transformation (RT) is a paradigmatic evolution of chronic lymphocytic leukemia (CLL) into a very aggressive large B-cell lymphoma conferring a dismal prognosis. The mechanisms driving RT remain largely unknown. We have characterized the whole genome, epigenome and transcriptome, combined with single-cell DNA/RNA sequencing analyses and functional experiments, of 19 CLL developing RT. Studying 54 longitudinal samples covering up to 19 years of disease course, we uncovered minute subclones carrying genomic, immunogenetic, and transcriptomic features of RT-cells already at CLL diagnosis, which were dormant for up to 19 years before transformation. We also identified new driver alterations, discovered a novel mutational signature (SBS-RT), recognized an oxidative phosphorylation (OXPHOS)high-B-cell receptor signaling (BCR)low transcriptional axis in RT, and showed that OXPHOS inhibition reduces the proliferation of RT cells. These findings demonstrate the early seeding of subclones driving advanced stages of cancer evolution and uncover therapeutic targets for the, once expanded, lethal RT.


This repository contains all the scripts, notebooks and reports to reproduce the scRNA-seq analysis of our paper "Early seeding of future Richter transformation in chronic lymphocytic leukemia", published in Nature Medicine in 2022. Here, we describe how to access the data, document the most important packages and versions used, and explain how to navigate the directories and files in this repository.


## Data

We have deposited the expression matrices, Seurat objects and metadata in Zenodo. You can download them following this link: [https://doi.org/10.5281/zenodo.6631966](https://doi.org/10.5281/zenodo.6631966), or executing the following commands in the
command line:

```{bash}
wget https://zenodo.org/record/6631966/files/Nadeu2022_scRNAseq.zip
unzip Nadeu2022_scRNAseq.zip
cd Nadeu2022_scRNAseq
```

The "Nadeu2022_scRNAseq" data directory comes with its own README.md, which give a full description of the dataset. To ease things, we also provide it here:


### Dataset description

The scRNA-seq data for this project was obtained for different patients (ids 12, 19, 63, 365, 3299) that underwent
chronic lymphoyctic leukemia to Richter transformation. For each patient we have sequential clinical samples: disease,
progression, treatment, Richter. The scRNA-seq data was sequenced in three different phases, which correspond to the
three subprojects one can find in the associated metadata:


#### Phase I: CLL_52

[Smart-seq2](https://www.nature.com/articles/nprot.2014.006) data for case 63.
This was a proof-of-concept we carried out early on for this project, to assess
if we could find Richter seeds in the initial time-points of the disease. Expression
matrices were called using [zUMIs](https://academic.oup.com/gigascience/article/7/6/giy059/5005022).


#### Phase II: BCLLATLAS_10

As the proof-of-concept worked well, we scaled the sequencing to all samples from
the other four cases: 12, 19, 365 and 3299. This corresponds to [10X Chromium v3](https://www.10xgenomics.com/products/single-cell-gene-expression) data, and was mapped with
[cellranger v4](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). For each case, we multiplexed the samples
coming from different disease timepoints with [cell hashing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1). This allowed us
to remove batch effects, detect doublets and decrease sequencing costs.


#### Phase III: BCLLATLAS_29

Because for some samples we obtained very few cells in the previous phase (low starting material),
we performed additionally 10X scRNA-seq for those samples: 365_07, 05_662, 00_30,
15-194, and 04-67.


### Metadata

#### Sequencing metadata

The file "richter_metadata.csv" has all the relevant information (subproject,
gem_id, library_id, library_name, type, and donor_id) for each sequencing library.
Sequencing libraries are derived from individual "GEM wells". As defined by
10X Genomcis: "GEM wells are a set of partitioned cells (Gel Beads-in-emulsion) from 
a single 10x Genomics Chromium™ Chip channel. One or more sequencing libraries can be
derived from a GEM well". More specifically, a GEM is: "an emulsion that contains a mixture
of biochemistry reagents (uniquely barcoded gel beads) and zero, one, or more
suspended cells/nuclei".

For cell-hashed libraries, each GEM well gives rise to two sequencing
libraries: cDNA (reverse-trascribed mRNAs) and hasthags oligonucleotides (HTO).
We give each GEM well a unique identifier, following the guideline described
[in this paper](https://academic.oup.com/gigascience/article/6/11/gix100/4557140?login=false).

All the information in this file helped us locate the right fastq files in our
high performance computer (HPC, cluster).


#### Feature barcode reference

For cell-hashed libraries, we need to know which sample is identified by each
barcode. This can be found in the "richter_feature_reference.csv" file,
and we used it as input for the "Feature Reference CSV", which needs to be
provided to cellranger as specificed in the [Antibody Capture analysis](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref)
vignette.


#### Sample metadata

The file "sample_metadata.csv" contains all the metadata for each sample. Each
sample represents a cell suspension obtained from a given disease time-point 
(i.e. diagnosis, progression, treatment, richter) for a given case (12, 19, 63,
365, 3299). The metadata includes donor_id, sample_id, tissue, time_point, 
treatment, etc.


#### CLL_52 metadata

See below.


### Expression matrices

For BCLLATLAS_10 and BCLLATLAS_29, we include the relevant files from the 
"outs" folder generated by cellranger:

* web_summary.html: html report with quality control metrics for the sequencing and mapping.
* metrics_summary.csv: csv file with the aforementioned QC metrics.
* filtered_feature_bc_matrix: a folder with 3 files, which correspond to the expresion matrix itself, the features and the cell barcodes.

For CLL_52, the expression matrix was called using [zUMIs](https://academic.oup.com/gigascience/article/7/6/giy059/5005022).
The output is a bit convoluted, so we have simplified it and save the expression
matrix as a csv file. In addition, the associated metadata can be found in the
"CLL_52_cell_metadata.csv".


### Seurat objects

In the Zenodo repository we have deposited all the initial, intermediate and final
[Seurat Objects](https://github.com/mojaveazure/seurat-object)
as .rds fies. Each object is the output of one of the analysis notebooks (.Rmd) that
one can find at the GitHub repository associated with this publication
(https://github.com/massonix/richter_transformation). We believe that this is
very helpful, as one can start exploring the data at any step of interest that
we show in our analysis.

In addition, we provide the output for the gene set enrichment analysis
(GSEA) that we obtained with [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/).
The file is called "patient_specific_gsea_raw_rt_vs_cll.rds" and before reading
it one should load clusterProfiler.


## Package versions

These are the versions of the most important packages used throughout all the analysis:

CRAN:

* [Seurat 4.0.3](https://satijalab.org/seurat/)
* dplyr 1.0.6
* purrr 0.3.4
* stringr 1.4.0
* tidyr 1.1.3
* ggplot2 3.3.3
* patchwork 1.1.1
* [harmony 1.0](https://github.com/immunogenomics/harmony)
* [UpSetR 1.4.0](https://github.com/hms-dbmi/UpSetR)


BioConductor:

* [clusterProfiler 3.18.1](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
* [infercnv 1.11.1](https://github.com/broadinstitute/infercnv)
* [GEOquery 2.62.1](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)


You can check the versions of other packages at the "Session Information" section of each html report.


## File system and name scheme

This repository contains 7 different analysis directories:

* 1-cellranger_mapping: scripts used to run cellranger in our cluster. It also contains QC metrics for different sequencing runs.
* 2-QC: quality control notebooks to assess the quality of sequencing and mapping, detect doublets, demultiplex hashtags, and filter poor-quality cells and genes.
* 3-clustering_and_annotation: we initially integrated all samples from all patients with [Harmony](), which allowed us to separate microenvironment from tumoral cells. Focusing on the latter, we then clustered and annotated all cells from each patient separately.
* 4-analysis_smart_seq2: notebooks with the full analysis of the Smart-seq2 data for patient 63.
* 5-clonal_evolution: notebooks that define the RT seed cells, and ensure that they are not artifacts.
* 6-differential_expression_analysis: differential expression analysis and gene set enrichment analysis between CLL and RT.
* 7-revision: notebooks to infer CNV and validate our findings with [an external dataset](https://aacrjournals.org/cancerdiscovery/article-abstract/11/12/3048/674669/Longitudinal-Single-Cell-Dynamics-of-Chromatin?redirectedFrom=fulltext).

The folder "bin" contains functions and variables that are used throughout different scripts. Finally, the folder "figures_and_tables" contains the scripts to generate the different figures derived from these analysis.


## Papers

If you want to learn more about Richter Syndrome I suggest reading the following reviews:

* [How we treat Richter syndrome](https://ashpublications.org/blood/article/123/11/1647/105732/How-we-treat-Richter-syndrome)
* [Biology and treatment of Richter syndrome](https://ashpublications.org/blood/article/131/25/2761/37138/Biology-and-treatment-of-Richter-syndrome)
