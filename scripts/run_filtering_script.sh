#!/bin/bash -l
#SBATCH --mem=80G
#SBATCH --array=6-9
#SBATCH -t 6-12:00:00
#SBATCH -p bigmemm

module load R/4.2

tmp_dir=$(mktemp -d /home/hemstrow/tempdir/filt.tmp.XXXXXXX)

Rscript snpR_filtering_script.R $SLURM_ARRAY_TASK_ID $tmp_dir
