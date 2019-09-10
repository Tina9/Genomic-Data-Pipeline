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

    mkdir -p indel_result
    mkdir -p indel_ori
    mkdir -p snp_result

    if((i%$number1==$number2))
    then

        $SAMTOOLS mpileup -B -q 2 -f $GENOME $normal_bam $tumor_bam | $VARSCAN somatic --output-snp ${tumor_sample}.snp.vcf --output-indel ${tumor_sample}.indel.vcf --mpileup 1 --min-coverage 10 --output-vcf 1
        $VARSCAN processSomatic ${tumor_sample}.indel.vcf

        grep -e "^#" -e "PASS" ${tumor_sample}.indel.Somatic.hc.vcf > indel_result/${tumor_sample}_indel_varscan.vcf

        mv ./${tumor_sample}.indel* indel_ori
        mv ./${tumor_sample}.snp.vcf snp_result
    fi
    i=$((i+1))

done
