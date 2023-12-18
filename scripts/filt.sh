#!/bin/bash -l
#SBATCH --array=43-46
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -A standby
#SBATCH --mem=10G
#SBATCH -J filter
#SBATCH --output=./slurm_files/slurm_%x_%A_%a.out
#SBATCH --error=./slurm_files/slurm_%x_%A_%a.err

module purge
module load r/4.2

infile=$1
outfile=$2
return_all_fst=$3
return_sfs=$4

parmfile=filt_parms.txt
ne_estimator_path=~/bin/Ne2-1L
par=4
tmpdir=/scratch/negishi/whemstro/tempdir

Rscript filt.R $SLURM_ARRAY_TASK_ID $infile $parmfile $outfile $return_all_fst $return_sfs $ne_estimator_path $par $tmpdir

