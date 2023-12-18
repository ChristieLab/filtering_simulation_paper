#!/bin/bash -l
#SBATCH --array=1-6
#SBATCH --mem=12G
#SBATCH -t 4:00:00
#SBATCH -J format_sweeps
#SBATCH -A standby

module load r/4.2

list=format_list.txt

string="sed -n ${SLURM_ARRAY_TASK_ID}p ${list}" 
str=$($string)

Rscript format_sweeps.R $str
