analysis_dir=$1
config_file=$2
number1=$3
number2=$4

ANNOVAR=/home/xzhang/softwares/Annovar/annovar/table_annovar.pl
HUMANDB=/home/xzhang/softwares/Annovar/annovar/humandb
STRELKA2ANNOVAR=/home/xzhang/Projects/Ping/5.snv_annovar/StreMut/strelka2annovar.py

cat $config_file | while read id
do
    arr=($id)
    sample=${arr[0]}
    snp_vcf=${arr[1]}
    indel_vcf=${arr[2]}

    if((i%$number1==$number2))
    then

    cp $indel_vcf ${sample}.indels.vcf.gz
    gunzip ${sample}.indels.vcf.gz
    grep -e "^#" -e "PASS" ${sample}.indels.vcf > ${sample}.filter.indels.vcf
    grep -v "#" ${sample}.filter.indels.vcf | awk '{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$10"\t"$11}' > ${sample}.indel.input

    $STRELKA2ANNOVAR ${sample}.indel.input ${sample}.indel.annovar.input

    $ANNOVAR $snp_vcf $HUMANDB -buildver hg38 -out ${sample}_snp -remove -protocol refGene,cytoBand,exac03,dbnsfp33a,avsnp150,cosmic70,1000g2015aug_all -operation g,r,f,f,f,f,f -nastring . -vcfinput
    $ANNOVAR ${sample}.indel.annovar.input $HUMANDB -buildver hg38 -out ${sample}_indel -remove -protocol refGene,cytoBand,exac03,dbnsfp33a,avsnp150,cosmic70,1000g2015aug_all -operation g,r,f,f,f,f,f -nastring .

    cut -f 1-88 ${sample}_snp.hg38_multianno.txt |awk -F "\t" '($6 == "exonic") && ($9 != "synonymous SNV") && ($9 != "unknown") {print $0}' > ${sample}_snp.filter.annovar.txt
    awk -F "\t" '($6 == "exonic") {print $0}' ${sample}_indel.hg38_multianno.txt > ${sample}_indel.filter.annovar.txt

    cat ${sample}_snp.filter.annovar.txt ${sample}_indel.filter.annovar.txt > ${sample}.filter.annovar.txt
    rm ${sample}.indels.vcf
    rm ${sample}.indel.input

    fi
    i=$((i+1))
done
