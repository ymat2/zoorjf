#!/bin/bash

#SBATCH --mem 16G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# raw_data1=~/raw_data/RJF_ABRC
raw_data2=~/raw_data/RJF_KMT
raw_data3=~/raw_data/RJF_iZoo
raw_data4=~/raw_data/RJF_TM
raw_data5=~/raw_data/RJF_MRYM
proj=~/RJF

cd ${proj}

for d in ${raw_data4} ${raw_data5}; do
  samples=($(ls ${d}))
  for sample in ${samples[@]}; do
    [ ! -e uncleaned/${sample} ] && mkdir -p uncleaned/${sample}
    paired_1=($(ls ${d}/${sample}/*_1.fq.gz))
    paired_2=($(ls ${d}/${sample}/*_2.fq.gz))
    cat ${paired_1[@]} > uncleaned/${sample}/${sample}_1.fq.gz
    cat ${paired_2[@]} > uncleaned/${sample}/${sample}_2.fq.gz
  done
done
