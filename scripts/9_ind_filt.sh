#!/bin/bash
#SBATCH --job-name=filter
#SBATCH -p high
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -t 4:00:00
#SBATCH --output=./slurm_filter_%A_%a.out
#SBATCH --error=./slurm_filter_%A_%a.err

module load vcftools

vcftools --gzvcf ../data/NUB_GBT_ypp.vcf.gz \
        --minGQ 13.0103 \
	--recode \
	--recode-INFO-all \
        --out ../data/NUB_GBT_yp_ggvcfs_hq





