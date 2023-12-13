samplemap=$1
bedfile=$2
mem=$3
tmp_dir=$4
par=$5
javapath=$6
gatkpath=$7
batchsize=$8

$javapath -jar -Xmx${mem}g -Xms${mem}g $gatkpath GenomicsDBImport \
       --genomicsdb-workspace-path ${bedfile}_db \
       --batch-size $batchsize \
       -L ${bedfile} \
       --sample-name-map ${samplemap} \
       --tmp-dir $tmp_dir \
       --reader-threads ${par}

