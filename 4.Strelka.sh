if [ $# != 4 ]; then
    echo  "USAGE: ./somatic.sh <analysis_dir> <config_file> <number1> <number2>"
   exit
fi

analysis_dir=$1
config_file=$2
number1=$3
number2=$4

########################################################
## point out the path of each software
GENOME=/home/xzhang/Databases/gatk_hg38/genome/Homo_sapiens_assembly38.fasta
INDEX=/home/xzhang/Databases/gatk_hg38/bwa_index/gatk_hg38
bed=/home/xzhang/Databases/gatk_hg38/annotation/CCDS/exon_probe.GRCh38.gene.150bp.bed
DBSNP=/home/xzhang/Databases/gatk_hg38/annotation/dbsnp_146.hg38.vcf.gz
kgSNP=/home/xzhang/Databases/gatk_hg38/annotation/1000G_phase1.snps.high_confidence.hg38.vcf.gz
kgINDEL=/home/xzhang/Databases/gatk_hg38/annotation/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

STRELKA=/home/xzhang/softwares/strelka-2.9.10.centos6_x86_64/bin
BCFTOOLS=/home/xzhang/anaconda3/bin/bcftools
# cd /home/haitaowang/Project/Genomics
# mkdir -p tmp
TMPDIR=/home/xzhang/Projects/Ping/4.Strelka/tmp

# mkdir -p alignment
mkdir -p $TMPDIR/Strelka_Result
cat $config_file |while read id
do
    arr=($id)
    normal_bam=${arr[1]}
    tumor_bam=${arr[2]}
    sample=${arr[0]}

    if((i%$number1==$number2))
    then

        if [ ! -f ${sample}_marked_fixed.bam ]; then

#####################################################
################ Step 1 : strelka2 ##################
#####################################################
cd $TMPDIR
mkdir -p ${sample}_somatic
$STRELKA/configureStrelkaSomaticWorkflow.py \
--normalBam $normal_bam \
--tumorBam $tumor_bam \
--referenceFasta $GENOME \
--runDir $TMPDIR/${sample}_somatic \
--exome \
--disableEVS \
--reportEVSFeatures \
--config=$STRELKA/configureStrelkaSomaticWorkflow.py.ini \
--snvScoringModelFile=/home/xzhang/softwares/strelka-2.9.10.centos6_x86_64/share/config/somaticSNVScoringModels.json \
--indelScoringModelFile=/home/xzhang/softwares/strelka-2.9.10.centos6_x86_64/share/config/somaticIndelScoringModels.json \
--outputCallableRegions

# execution on a single local machine with 20 parallel jobs
$TMPDIR/${sample}_somatic/runWorkflow.py -m local -j 4

################################################################
################ Step 2: merge snv and indel result#############
################################################################

$BCFTOOLS concat -a $TMPDIR/${sample}_somatic/results/variants/somatic.snvs.vcf.gz $TMPDIR/${sample}_somatic/results/variants/somatic.indels.vcf.gz -o $TMPDIR/Strelka_Result/${sample}_strelka.vcf


grep -e "^#" -e "PASS" $TMPDIR/Strelka_Result/${sample}_strelka.vcf > $TMPDIR/Strelka_Result/${sample}_strelka_filter.vcf

        fi
    fi
    i=$((i+1))
done
