if [ $# != 4 ]; then
    echo  "USAGE: ./control-freec <analysis_dir> <config_file> <number1> <number2>"
    exit
fi

analysis_dir=$1
config_file=$2
number1=$3
number2=$4

BEDTOOLS="/home/xzhang/anaconda3/bin/bedtools"
REF_CNV="/home/xzhang/Projects/Ping/6.FREEC/6.3.cnv_total/ref/chromosome.hg38.sorted.txt"
cat $config_file |while read id
do
    arr=($id)
    sample=${arr[0]}
    cnv=${arr[1]}

    if((i%$number1==$number2))
    then
    awk '$4 != -1 {print $0}' $cnv | cut -f 1,4,6 | awk -F '[\t:-]' '{print "chr"$3"\t"$4"\t"$5"\t"$2}' | grep -v chrGene | sort -k1,1 -k2,2n > ${sample}_cnv.bed
    $BEDTOOLS intersect -wa -wb -a $REF_CNV -b ${sample}_cnv.bed -loj| awk '{print $1"-"$2"-"$3"\t"$9}' | awk '{sum2[$1] += $2; count2[$1]++}; END{for (id in sum2) {print id, sum2[id]/count2[id]}}' | sed 's/ /\t/g' > ${sample}_cnv.txt

    fi
    i=$((i+1))

done
