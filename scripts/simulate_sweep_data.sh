#!/bin/bash -l
#SBATCH --mem=34G
#SBATCH -t 4-4:00:00
#SBATCH -J sim
#SBATCH -A christ99
#SBATCH --array=79-114

module load r/4.2

parmfile=sweep_sim_parms.txt

Rscript simulate_sweep_data.R $parmfile $SLURM_ARRAY_TASK_ID /home/whemstro/bin/msms3.2rc-b163.jar /home/whemstro/bin/jre1.8.0_371/bin/java
