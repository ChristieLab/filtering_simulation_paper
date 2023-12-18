#!/bin/bash -l
#SBATCH --mem=100G
#SBATCH -t 5-8:00:00
#SBATCH -J format
#SBATCH -A christ99

module load r/4.2

Rscript format.R humans_merged_snps_poly.vcf.gz sample_meta.txt ../snpR_RDS/humans.RDS
