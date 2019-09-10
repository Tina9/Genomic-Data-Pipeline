if [ $# != 4 ]; then
    echo  "USAGE: ./control-freec <analysis_dir> <config_file> <number1> <number2>"
    exit
fi

analysis_dir=$1
config_file=$2
number1=$3
number2=$4

PYTHON=~/anaconda3/bin/python

cat $config_file |while read id
do
    arr=($id)
    normal_bam=${arr[1]}
    tumor_bam=${arr[2]}
    sample=${arr[0]}

    if((i%$number1==$number2))
    then

    start=$(date +%s.%N)
    echo "Control-Freec" `date`
    $PYTHON proconf.py $tumor_bam $normal_bam $sample
    echo "Control-Freec" `date`
    dur=$(echo "$(date +%s.%N) - $start" | bc)
    printf "Execution time for Control-Freec : %.6f seconds" $dur
    echo

    fi
    i=$((i+1))

done
