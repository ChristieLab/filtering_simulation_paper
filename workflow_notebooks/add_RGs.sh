bamfile=$1 # bamfile, without the final extension (.bam)
fastq=$2 # one of the fastq files from which the bam file was produced
basename=$3 # the basename (base bam name) for the bamfile, without any extensions
javapath=$4
picardpath=$5
platform=$6

# find the PU info
fastq_head="$(head -n 1 $fastq)"
IFS=':' read -ra rinf <<< "$fastq_head"

$javapath -jar $picardpath AddOrReplaceReadGroups \
        I=${bamfile}.bam \
        O=${bamfile}.RG.bam \
        RGLB=Lib-${basename} \
        RGPL=${platform} \
        RGPU=${platform}_${basename} \
        RGSM=${basename} \
        VALIDATION_STRINGENCY=LENIENT

samtools index ${bamfile}.RG.bam
