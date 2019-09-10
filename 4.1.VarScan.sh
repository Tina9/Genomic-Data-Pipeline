if [ $# != 4 ]; then
    echo  "USAGE: ./varscan.sh <analysis_dir> <config_file> <number1> <number2>"
    exit
fi

analysis_dir=$1
config_file=$2
number1=$3
number2=$4

GENOME=/home/xzhang/Databases/gatk_hg38/genome/Homo_sapiens_assembly38.fasta
VARSCAN=/home/xzhang/anaconda3/envs/VarScan/bin/varscan
SAMTOOLS=/home/xzhang/anaconda3/bin/samtools
BGZIP=/home/xzhang/anaconda3/bin/bgzip
TABIX=/home/xzhang/anaconda3/bin/tabix
BCFTOOLS=/home/xzhang/anaconda3/bin/bcftools

source activate VarScan

cat $config_file | while read id
do
    arr=($id)
    normal_bam=${arr[1]}
    tumor_bam=${arr[2]}
    normal_sample=$(basename $normal_bam | cut -d"_" -f 1)
    tumor_sample=${arr[0]}

    if((i%$number1==$number2))
    then
        if [ ! -f "${normal_sample}.mpileup" ]; then

            $SAMTOOLS mpileup -f $GENOME $normal_bam > ${normal_sample}.mpileup

        fi
            $SAMTOOLS mpileup -f $GENOME $tumor_bam > ${tumor_sample}.mpileup
            $VARSCAN somatic ${normal_sample}.mpileup ${tumor_sample}.mpileup $tumor_sample --output-vcf 1 --min-coverage 10
            $VARSCAN

            rm ${normal_sample}.mpileup
            rm ${tumor_sample}.mpileup

            $BGZIP ${tumor_sample}.snp.vcf
            $TABIX ${tumor_sample}.snp.vcf.gz

            $BGZIP ${tumor_sample}.indel.vcf
            $TABIX ${tumor_sample}.indel.vcf.gz

            $BCFTOOLS concat -a ${tumor_sample}.snp.vcf.gz ${tumor_sample}.indel.vcf.gz -o ${tumor_sample}_varscan.vcf

            grep -e "^#" -e "PASS" ${tumor_sample}_varscan.vcf > ${tumor_sample}_varscan_filter.vcf

    fi
    i=$((i+1))

done

conda deactivate
