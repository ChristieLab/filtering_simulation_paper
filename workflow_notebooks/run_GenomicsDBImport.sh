samplemap=$1
bedfile=$2
mem=$3
tmp_dir=$4
par=$5

basebed=$(basename ${bedfile})

gatk --java-options "-Xmx${mem}g -Xms${mem}g" GenomicsDBImport  \
       --genomicsdb-workspace-path ${bedbed}_db \
       --batch-size 15 \
       -L ${bedfile} \
       --sample-name-map ${samplemap} \
       --tmp-dir $tmp_dir \
       --reader-threads ${par}


echo "build database done at: `date`"
echo ""

