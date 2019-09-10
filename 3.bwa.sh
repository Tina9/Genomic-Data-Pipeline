
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

source /etc/profile
source /etc/profile.d/modules.sh
module add singularity/2.6.1
export SIMG=/share/apps/singularity/simg/bioconda/bioconda-anaconda3-gwlab-2018.12
SINGULARITY_EXEC="singularity exec --bind /data:/data,/scratch:/scratch,/share:/share,/home:/home ${SIMG} bash -c "

home="/home/xzhang"
scratchHome="/scratch/xzhang"

GENOME=/home/xzhang/Databases/gatk_hg38/genome/Homo_sapiens_assembly38.fasta
INDEX=/home/xzhang/Databases/gatk_hg38/bwa_index/gatk_hg38
bed=/home/xzhang/Databases/gatk_hg38/annotation/CCDS/exon_probe.GRCh38.gene.150bp.bed
DBSNP=/home/xzhang/Databases/gatk_hg38/annotation/dbsnp_146.hg38.vcf.gz
kgSNP=/home/xzhang/Databases/gatk_hg38/annotation/1000G_phase1.snps.high_confidence.hg38.vcf.gz
kgINDEL=/home/xzhang/Databases/gatk_hg38/annotation/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
GATK=/opt/conda/bin/gatk
BWA=/opt/conda/bin/bwa
SAMTOOLS=/opt/conda/envs/samtools/bin/samtools
###############################
cat $config_file | while read id
do
    arr=($id)
    fq1=${arr[1]}
    fq2=${arr[2]}
    sample=${arr[0]}

    if((i%$number1==$number2))
    then

        if [  ! -f ${sample}_recal.bam ]; then
#####################################################
################ Step 1 : Alignment #################
#####################################################
start=$(date +%s.%N)
echo $BWA `date`
${SINGULARITY_EXEC}"$BWA mem -t 4 -M -R \"@RG\tID:$sample\tSM:$sample\tLB:WES\tPL:Illumina\" $INDEX $fq1 $fq2 1>$sample.sam"
echo $BWA `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for BWA : %.6f seconds" $dur
echo


####################################################
############### Step 2: Sort and Index #############
####################################################
start=$(date +%s.%N)
echo SortSam `date`
${SINGULARITY_EXEC}"$GATK --java-options \"-Xmx25G -Djava.io.tmpdir=./\"  SortSam \
-SO coordinate -I $sample.sam -O $sample.bam #1>log.sort 2>&1"
${SINGULARITY_EXEC}"source activate samtools; $SAMTOOLS index $sample.bam"
echo SortSam `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for SortSam : %.6f seconds" $dur
echo
rm $sample.sam


#####################################################
################ Step 3: Basic Statistics ###########
#####################################################
start=$(date +%s.%N)
echo stats `date`
${SINGULARITY_EXEC}"source activate samtools; $SAMTOOLS flagstat $sample.bam > ${sample}.alignment.flagstat"
${SINGULARITY_EXEC}"source activate samtools; $SAMTOOLS stats $sample.bam > ${sample}.alignment.stat"
echo plot-bamstats -p ${sample}_QC  ${sample}.alignment.stat
echo stats `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for Basic Statistics : %.6f seconds" $dur
echo


####################################################
###### Step 4: multiple filtering for bam files ####
####################################################

###MarkDuplicates###
start=$(date +%s.%N)
echo MarkDuplicates `date`
${SINGULARITY_EXEC}"$GATK  --java-options \"-Xmx25G -Djava.io.tmpdir=./\" MarkDuplicates  \
-I $sample.bam -O ${sample}_marked.bam -M $sample.metrics --REMOVE_DUPLICATES true"
echo MarkDuplicates `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for MarkDuplicates : %.6f seconds" $dur
echo
rm $sample.bam
rm $sample.bam.bai


###FixMateInfo###
start=$(date +%s.%N)
echo FixMateInfo `date`
${SINGULARITY_EXEC}"$GATK --java-options \"-Xmx25G -Djava.io.tmpdir=./\" FixMateInformation \
-I ${sample}_marked.bam -O ${sample}_marked_fixed.bam -SO coordinate"
${SINGULARITY_EXEC}"source activate samtools; $SAMTOOLS index ${sample}_marked_fixed.bam"
echo FixMateInfo `date`
dur=$(echo "$(date +%s.%N) - $start" | bc)
printf "Execution time for FixMateInfo  : %.6f seconds" $dur
echo
rm ${sample}_marked.bam

####################################################
################### Step 5: recal ##################
####################################################
bam=${sample}_marked_fixed.bam
${SINGULARITY_EXEC}"$GATK  --java-options \"-Xmx20G -Djava.io.tmpdir=./\"   BaseRecalibrator \
-I $bam -R $GENOME --output ${sample}_recal.table --known-sites $kgSNP --known-sites $kgINDEL"
${SINGULARITY_EXEC}"$GATK  --java-options \"-Xmx20G -Djava.io.tmpdir=./\"   ApplyBQSR \
-I $bam -R $GENOME --output ${sample}_recal.bam -bqsr ${sample}_recal.table"
${SINGULARITY_EXEC}"source activate samtools; $SAMTOOLS index ${sample}_recal.bam"
rm ${sample}_marked_fixed.bam
rm ${sample}_marked_fixed.bam.bai
rm ${sample}_recal.bai

        fi
    fi
    i=$((i+1))
done
