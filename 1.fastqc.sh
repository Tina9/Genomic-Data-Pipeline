#### Add tasks submission commands ###

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

fastqc="/opt/conda/bin/fastqc"


cat $config_file | while read id
do
    arr=($id)
    fq1=${arr[1]}
    fq2=${arr[2]}
    sample=${arr[0]}

    if((i%$number1==$number2));
    then
    ${SINGULARITY_EXEC}"$fastqc -t 2 $fq1 $fq2"
    fi
    i=$((i+1))
done
