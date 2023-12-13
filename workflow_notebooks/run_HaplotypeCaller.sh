bamfile=$1
reference=$2
tmp_dir=$3
mem=$4
javapath=$5
gatkpath=$6

$javapath -jar -Xmx${mem}g -Xms${mem}g -Djava.io.tmpdir=${tmp_dir} $gatkpath HaplotypeCaller \
        -ERC GVCF \
        -R $reference \
        -I $bamfile \
        -O $bamfile.hapcalls.gvcf.gz
