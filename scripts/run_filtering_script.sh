#!/bin/bash -l
#SBATCH --mem=50G
#SBATCH --array=1-14
#SBATCH -t 6-12:00:00
#SBATCH -p med

module load R/4.2

tmp_dir=$(mktemp -d /home/hemstrow/tempdir/filt.tmp.XXXXXXX)

Rscript snpR_filtering_script.R $SLURM_ARRAY_TASK_ID $tmp_dir
