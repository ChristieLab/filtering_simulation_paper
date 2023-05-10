#!/bin/bash -l
#SBATCH --mem=80G
#SBATCH --array=27-35
#SBATCH -t 6-12:00:00
#SBATCH -p bigmemm

module load R/4.2

parmfile=monarch_flt_parmfile.txt
outprefix=../results/filtering_results_monarchs
ne_estimator_path=/home/hemstrow/bin/Ne2-1L

tmp_dir=$(mktemp -d /home/hemstrow/tempdir/filt.tmp.XXXXXXX)

Rscript snpR_filtering_script.R \
  $SLURM_ARRAY_TASK_ID \
  $tmp_dir \
  $parmfile \
  ${outprefix}_${SLURM_ARRAY_TASK_ID}.RDS \
  $ne_esttimator_path