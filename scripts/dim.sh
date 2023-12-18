#!/bin/bash -l
#SBATCH --array=1-11
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -A standby
#SBATCH --mem=220G
#SBATCH -J filter

module purge
module load r/4.2

list=file_list.txt
string="sed -n ${SLURM_ARRAY_TASK_ID}p ${list}"
str=$($string)

Rscript dim.R $str
