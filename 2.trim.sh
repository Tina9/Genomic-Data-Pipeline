### There is a great need to do fastqc after trim

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

########################################################
## point out the path of each software
bin_trim_galore="/opt/conda/envs/py36/bin/trim_galore"

mkdir -p clean_fastq

cat $config_file |while read id
do
    arr=($id)
    fq1=${arr[1]}
    fq2=${arr[2]}
    sample=${arr[0]}

    if((i%$number1==$number2))
    then
        if [ ! -f $analysis_dir/clean_fastq/${fq2_base}_val_2.fq.gz ]; then
        echo "start trim_galore for $sample" `date`
        ${SINGULARITY_EXEC}"source activate py36; $bin_trim_galore -q 25 --phred33 --length 50 -e 0.1 --stringency 3 --paired -o $analysis_dir/clean_fastq $fq1 $fq2"
        echo "end trim_galore for $sample" `date`
        fi
    fi
    i=$((i+1))
done
########################################################
