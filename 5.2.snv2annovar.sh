if [ $# != 4 ]; then
    echo  "USAGE: ./trim.sh <analysis_dir> <config_file> <number1> <number2>"
    exit
fi

#############################################################
#############################################################
analysis_dir=$1
config_file=$2
number1=$3
number2=$4

ANNOVAR=/home/xzhang/softwares/Annovar/annovar/table_annovar.pl
HUMANDB=/home/xzhang/softwares/Annovar/annovar/humandb
TABIX=/home/xzhang/anaconda3/bin/tabix
BCFTOOLS=/home/xzhang/anaconda3/bin/bcftools
BGZIP=/home/xzhang/anaconda3/bin/bgzip

cat $config_file | while read id
do
    arr=($id)
    sample=${arr[0]}
    snp_vcf=${arr[1]}
    indel_vcf=${arr[2]}

    if((i%$number1==$number2))
    then

    $ANNOVAR $snp_vcf $HUMANDB -buildver hg38 -out ${sample}_snp -remove -protocol refGene,cytoBand,exac03,dbnsfp33a,avsnp150,cosmic70,1000g2015aug_all -operation g,r,f,f,f,f,f -nastring . -vcfinput

    cut -f 1-88 ${sample}_snp.hg38_multianno.txt |awk -F "\t" '($6 == "exonic") && ($9 != "synonymous SNV") && ($9 != "unknown") {print $0}' > ${sample}_snp.filter.annovar.txt

    cat ${sample}_snp.filter.annovar.txt ${sample}_indel.filter.annovar.txt > ${sample}.filter.annovar.txt

    fi
    i=$((i+1))
done
