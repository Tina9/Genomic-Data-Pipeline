if [ $# != 4 ]; then
    echo  "USAGE: ./control-freec <analysis_dir> <config_file> <number1> <number2>"
    exit
fi

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
CNV_PLOT="/home/xzhang/Projects/Ping/6.FREEC/6.2.merge/single_CNV.R"

cat $config_file |while read id
do
    arr=($id)
    tumor_cnv=${arr[1]}
    organ_cnv=${arr[2]}
    sample=${arr[0]}

    if((i%$number1==$number2))
    then
    grep '' $tumor_cnv $organ_cnv | grep -v Chromosome | awk '$4 != -1 {print $0}'  | cut -d "/" -f 9 | sed s/_recal.bam_ratio.txt:/\\t/g > ${sample}_CN.txt
    ${SINGULARITY_EXEC}"Rscript $CNV_PLOT ${sample}_CN.txt $sample"

    fi
    i=$((i+1))

done
