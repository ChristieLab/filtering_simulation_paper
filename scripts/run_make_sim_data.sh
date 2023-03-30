#!/bin/bash -l
#SBATCH --mem=20G
#SBATCH -n 30
#SBATCH -N 1
#SBATCH -t 7-12:00:00
#SBATCH -J make_sim
#SBATCH -A beagle
#SBATCH --output=slurm_make_sim_%A_%a.out
#SBATCH --error=slurm_make_sim_%A_%a.err

module load r/4.1

Rscript /scratch/bell/whemstro/filtering_simulation_paper/scripts/make_sim_data.R
