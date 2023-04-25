#!/bin/bash -l
#SBATCH --mem=80G
#SBATCH --array=1-14
#SBATCH -t 6-12:00:00
#SBATCH -A beagle

module load r/4.2

tmp_dir=$(mktemp -d /scratch/bell/whemstro/tempdir/filt.tmp.XXXXXXX)

Rscript snpR_filtering_script_ypp.R $SLURM_ARRAY_TASK_ID $tmp_dir
