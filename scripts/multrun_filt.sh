#!/bin/bash -l
#SBATCH -t 2:00:00
#SBATCH --array=1-20
#SBATCH --mem=1G
#SBATCH -J filter_handler
#SBATCH -A standby
#SBATCH --output=./slurm_files/slurm_%x_%A_%a.out
#SBATCH --error=./slurm_files/slurm_%x_%A_%a.err

set -e

joblist=filt_jobs.txt

echo "Starting!"

string="sed -n ${SLURM_ARRAY_TASK_ID}p ${joblist}"
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2, $3, $4, $5}')
set -- $var
infile=$1
outfile=$2
memory=$3
return_all_fst=$4
return_sfs=$5

echo "Running job on: ${infile}"
echo "with ${memory}G of memory"
echo "to outfile: ${outfile}"
echo "Saving per-locus fst: ${return_all_fst}"
echo "Calculating SFS: ${return_sfs}"

# update the required memory in the job
cp filt.sh filt_r${SLURM_ARRAY_TASK_ID}.sh
perl -i -pe "s/#SBATCH --mem=.+/#SBATCH --mem=${memory}G/g" filt_r${SLURM_ARRAY_TASK_ID}.sh

# run the jobs
j=$(sbatch --parsable filt_r${SLURM_ARRAY_TASK_ID}.sh $infile $outfile $return_all_fst $return_sfs)

echo "Job ID: ${j}"

head -n 10 filt_r${SLURM_ARRAY_TASK_ID}.sh

rm filt_r${SLURM_ARRAY_TASK_ID}.sh
