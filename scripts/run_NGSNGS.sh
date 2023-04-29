#!/bin/bash -l
#SBATCH --mem=40G
#SBATCH -t 6-12:00:00
#SBATCH -A mcclintock
#SBATCH -n 64
#SBATCH -N 1
#SBATCH --array=1-1

ref=
variants=
reads=1000
seed=6216
read_legnth=150
output_prefix=test_reads

~/bin/NGSNGS/ngsngs -i $ref -r $reads -t 64 -s $seed -cl $read_length \
	-seq PE -f fq.gz \
	-o $output_prefix \
	-vcf $variants \
	-id ${SLURM_ARRAY_TASK_ID}

