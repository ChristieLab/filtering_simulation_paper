vcf=$1
ref=$2
QD="${3}"
FS="${4}"
SOR="${5}"
MQ="${6}"
MQRS="${7}"
RPRS="${8}"
GQ="${9}"
mem=${10}
javapath=${11}
gatkpath=${12}

$javapath -jar -Xmx${mem}g -Xms${mem}g $gatkpath VariantFiltration \
       -R $ref \
       -V ${vcf}.vcf \
       --filter-name "QDf" \
       --filter-expression "QD < $QD" \
       --filter-name "FSf" \
       --filter-expression "FS > $FS" \
       --filter-name "SORf" \
       --filter-expression "SOR > $SOR" \
       --filter-name "MQf" \
       --filter-expression "MQ < $MQ" \
       --filter-name "MQRSf" \
       --filter-expression "MQRankSum < $MQRS" \
       --filter-name "RPRSf" \
       --filter-expression "ReadPosRankSum < $RPRS" \
       -O hard_filt_${vcf}.vcf

bcftools view --types snps -m 2 -M 2 hard_filt_${vcf}.vcf > hard_filt_temp_${vcf}.vcf

vcftools --vcf hard_filt_temp_${vcf}.vcf \
        --remove-filtered-all \
        --remove-indels \
        --minGQ $GQ \
        --recode \
        --recode-INFO-all \
        --out hard_filt_pass_${vcf}

rm hard_filt_temp_${vcf}.vcf
