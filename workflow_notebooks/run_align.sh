c1=$1 # A file containing the names of the files containing data for each individual (RA). The .fastq at the end is cutt off.
c2=$2 # A file containing the names of the files containing data for each individual (RA). The .fastq at the end is cutt off.
ref=$3 # the path to the genome for the alignment
mapQ=$4 # mapping quality min removing poorly mapped reads
rmdup=$5 # 0 or 1. Should suspected PCR duplicates be removed?
improperpairs=$6 # 0 or 1. Should improperly paired reads be removed?


bwa mem $ref ${c1}.fastq ${c2}.fastq | samtools view -Sb - | samtools sort - -n -o ${c1}.sort.bam # align and sort by name
samtools fixmate -r -m ${c1}.sort.bam ${c1}.fixmate.bam # fixmate, marks some info
samtools sort -o ${c1}.psort.bam ${c1}.fixmate.bam # sort by position

# filter PCR dups if requested, then do map quality filtering
if [[ $rmdup -eq 1 ]]
then
  samtools markdup -r ${c1}.psort.bam ${c1}.markdup.bam # remove dups, requires fixmate
  samtools view -q $mapQ -b ${c1}.markdup.bam > ${c1}.q1.bam # remove poorly mapped
  # rm ${c1}.markdup.bam
else
  samtools view -q $mapQ -b ${c1}.psort.bam > ${c1}.q1.bam # remove poorly mapped
fi

# resort and refixmate
samtools sort -n -o ${c1}.namesort.bam ${c1}.q1.bam # sort by name again
samtools fixmate -m ${c1}.namesort.bam ${c1}.fixmate.bam # mark bad mates again

# improper pairs filtering if requested, then resort
if [[ $improperpair -eq 1 ]]
then
  samtools view -f 0x2 -b ${c1}.fixmate.bam > ${c1}.flt.bam # remove improper pairs
  samtools sort -o ${c1}.sort.flt.bam ${c1}.flt.bam # sort
else
  samtools sort -o ${c1}.sort.flt.bam ${c1}.fixmate.bam # sort
fi

# index
samtools index ${c1}.sort.flt.bam

# clean
if [[ $rmdup -eq 1 ]]
then
  rm ${c1}.markdup.bam
fi

if [[ $improperpair -eq 1 ]]
then
  ${c1}.flt.bam
fi

rm ${c1}.fixmate.bam
rm ${c1}.q1.bam
rm ${c1}.psort.bam
rm ${c1}.namesort.bam
