#!/bin/bash -l
#SBATCH --mem=100G
#SBATCH --array=1-35
#SBATCH -t 6-12:00:00
#SBATCH -p bigmemm

module load R/4.2

parmfile=ypp_filt_parmfile.txt
outprefix=../results/filtering_results_ypp
ne_estimator_path=/home/hemstrow/bin/Ne2-1L

tmp_dir=$(mktemp -d /home/hemstrow/tempdir/filt.tmp.XXXXXXX)

Rscript snpR_filtering_script_ypp.R \
  $SLURM_ARRAY_TASK_ID \
  $tmp_dir \
  $parmfile \
  ${outprefix}_${SLURM_ARRAY_TASK_ID}.RDS \
  $ne_estimator_path
