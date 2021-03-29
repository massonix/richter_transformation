#!/bin/bash
    
#SBATCH --job-name="yiu7os1e_xng0z3m5"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/yiu7os1e_xng0z3m5_%x_%J.err
#SBATCH --output=./log/yiu7os1e_xng0z3m5_%x_%J.out
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/projects/BCLLATLAS_29/jobs/yiu7os1e_xng0z3m5/fastq --id yiu7os1e_xng0z3m5 --chemistry SC3Pv3 --expect-cells 5000 --localcores 24 --localmem 64 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;


echo [`date "+%Y-%m-%d %T"`] job finished