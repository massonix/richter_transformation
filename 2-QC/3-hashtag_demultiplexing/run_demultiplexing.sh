for i in $(grep hashed_cdna ../../1-cellranger_mapping/data/richter_metadata.csv);
do
	subproject=$(echo $i | cut -d',' -f1);
	gem_id=$(echo $i | cut -d',' -f2);
	library_name=$(echo $i | cut -d',' -f4);
	echo $subproject;
	echo $gem_id;
	echo $library_name;
	R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc'); \
	library(BiocStyle); \
	path_to_knit <- '/home/rmassonix/Desktop/PhD/RICHTER/'; \
	subproject <- '${subproject}'; \
	gem_id <- '${gem_id}'; \
	library_name <- '${library_name}'; \
	rmd <- ifelse(library_name %in% c('ICGC_3299_01', 'ICGC_3299_02'), '01-demultiplex_hashtags_ICGC_3299.Rmd', '01-demultiplex_hashtags.Rmd'); \
	save_object_path <- 'results/R_objects/demultiplexed/seurat_${gem_id}_demultiplexed.rds'; \
	rmarkdown::render( \
		rmd, \
		output_file = paste('reports/demultiplex_hashtags-', gem_id, '.html', sep = ''), \
		knit_root_dir = path_to_knit \
	)";

done
