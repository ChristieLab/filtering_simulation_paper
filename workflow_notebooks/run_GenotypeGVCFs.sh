bedfile=$1
ref=$2
mem=$3
tmp_dir=$4
javapath=$5
gatkpath=$6

# run the genotyping
echo "genotype start at: `date`"
$javapath -jar -Xmx${mem}g -Xms${mem}g $gatkpath GenotypeGVCFs \
  -R $ref \
  -L $bedfile \
  -V gendb://${bedfile}_db \
  -O raw_${bedfile}.vcf \
  --tmp-dir $tmp_dir \
  --new-qual \
  --max-alternate-alleles 2 \
  --disable-bam-index-caching
