#!/bin/bash
    
#SBATCH --job-name="of7uaqm6_zbhy4tn9"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/of7uaqm6_zbhy4tn9_%x_%J.err
#SBATCH --output=./log/of7uaqm6_zbhy4tn9_%x_%J.out
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --fastqs /scratch/devel/rmassoni/RICHTER/current/1-cellranger_mapping/projects/BCLLATLAS_29/jobs/of7uaqm6_zbhy4tn9/fastq --id of7uaqm6_zbhy4tn9 --chemistry SC3Pv3 --expect-cells 5000 --localcores 24 --localmem 64 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;


echo [`date "+%Y-%m-%d %T"`] job finished