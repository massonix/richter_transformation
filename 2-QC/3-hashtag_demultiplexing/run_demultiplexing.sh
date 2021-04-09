for i in $(grep hashed_cdna ../../1-cellranger_mapping/data/richter_metadata.csv);
do
	subproject=$(echo $i | cut -d',' -f1);
	gem_id=$(echo $i | cut -d',' -f2);
	echo $subproject;
	echo $gem_id;
	R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc'); \
	library(BiocStyle); \
	path_to_knit <- '/home/rmassonix/Desktop/PhD/RICHTER/'; \
	subproject <- '${subproject}'; \
	gem_id <- '${gem_id}'; \
	save_object_path <- 'results/R_objects/demultiplexed/seurat_${gem_id}_demultiplexed.rds'; \
	rmarkdown::render( \
		'01-demultiplex_hashtags.Rmd', \
		output_file = paste('reports/demultiplex_hashtags-', gem_id, '.html', sep = ''), \
		knit_root_dir = path_to_knit \
	)";

done
