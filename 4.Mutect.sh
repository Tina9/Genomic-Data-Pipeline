analysis_dir=$1
config_file=$2
number1=$3
number2=$4

source /etc/profile
source /etc/profile.d/modules.sh
module add singularity/2.6.1
export SIMG=/share/apps/singularity/simg/bioconda/mutect-1.1.7.simg
SINGULARITY_EXEC="singularity exec --bind /data:/data,/scratch:/scratch,/share:/share,/home:/home ${SIMG} bash -c "

home="/home/xzhang"
scratchHome="/scratch/xzhang"

GENOME=/home/xzhang/Databases/gatk_hg38/genome/Homo_sapiens_assembly38.fasta


cat $config_file | while read id
do
    arr=($id)
    normal_bam=${arr[1]}
    tumor_bam=${arr[2]}
    sample=${arr[0]}

    if((i%$number1==$number2))
    then

    ${SINGULARITY_EXEC}"java -Xmx10g -jar /bin/mutect-1.1.7 \
    --analysis_type MuTect \
    --reference_sequence $GENOME \
    --input_file:normal $normal_bam \
    --input_file:tumor $tumor_bam \
    --vcf ${sample}_Mutect.vcf"

    grep -e "^#" -e "PASS" ${sample}_Mutect.vcf > ${sample}_Mutect_filter.vcf

    fi
    i=$((i+1))

done
