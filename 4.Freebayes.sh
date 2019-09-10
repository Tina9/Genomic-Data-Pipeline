if [ $# != 4 ]; then
    echo  "USAGE: ./trim.sh <analysis_dir> <config_file> <number1> <number2>"
    exit
fi
#  "USAGE: ./trim.sh <analysis_dir> <config_file> <number1> <number2>"
#  example : sbatch trim.sh ./ config 3 0

analysis_dir=$1
config_file=$2
number1=$3
number2=$4

FREEBAYES=/home/xzhang/anaconda3/bin/freebayes
GENOMIC=/home/xzhang/Databases/gatk_hg38/genome/Homo_sapiens_assembly38.fasta
TABIX=/home/xzhang/anaconda3/bin/tabix
BCFTOOLS=/home/xzhang/anaconda3/bin/bcftools
BGZIP=/home/xzhang/anaconda3/bin/bgzip

cat $config_file | while read id
do
    arr=($id)
    normal_bam=${arr[1]}
    tumor_bam=${arr[2]}
    normal_sample=$(basename $normal_bam | cut -d"_" -f 1)
    tumor_sample=${arr[0]}

    if((i%$number1==$number2))
    then

        if [ ! -f "${normal_sample}.vcf.gz" ]; then

            $FREEBAYES -C 5 -f $GENOMIC $normal_bam > ${normal_sample}.vcf
                        $BGZIP ${normal_sample}.vcf
            $TABIX ${normal_sample}.vcf.gz

        fi

        $FREEBAYES -C 5 -f $GENOMIC $tumor_bam > ${tumor_sample}.vcf

        $BGZIP ${tumor_sample}.vcf
        $TABIX ${tumor_sample}.vcf.gz

        $BCFTOOLS isec -C -p $tumor_sample ${tumor_sample}.vcf.gz ${normal_sample}.vcf.gz
        mv $tumor_sample/0000.vcf ${tumor_sample}_freebayes.vcf
        rm -r $tumor_sample
        rm ${tumor_sample}.vcf.gz ${tumor_sample}.vcf.gz.tbi
    fi
    i=$((i+1))
done
