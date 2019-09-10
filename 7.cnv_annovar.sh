if [ $# != 4 ]; then
    echo  "USAGE: ./7.annovar.sh <analysis_dir> <config_file> <number1> <number2>"
    exit
fi

analysis_dir=$1
config_file=$2
number1=$3
number2=$4

ANNOVAR=/home/xzhang/softwares/Annovar/annovar/table_annovar.pl
HUMANDB=/home/xzhang/softwares/Annovar/annovar/humandb

cat $config_file | while read id
do
    arr=($id)
    cnv_file=${arr[1]}
    sample=${arr[0]}

    if((i%$number1==$number2))
    then

    grep -v "Chromosome" $cnv_file | awk '($5 != -1 && $5 != 2 && $3 != -1) {print $0}' | awk -F '[:-]' '{print $1"\t"$2"\t"$3}' | awk '{print $6"\t"$7"\t"$8"\t"0"\t"0"\t"$4"\t"$5}' > ${sample}.csv

    $ANNOVAR ${sample}.csv $HUMANDB -buildver hg38 -out $sample -remove -protocol refGene,cytoBand -operation g,r -nastring NA -otherinfo

    #rm ${sample}.csv

    fi
    i=$((i+1))

done
