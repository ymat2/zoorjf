#!/bin/bash

#SBATCH -a 1-9
#SBATCH --mem 32G

proj=~/RJF
raw_data=${proj}/uncleaned
clean_data=${proj}/clean
samples=($(ls ${raw_data}))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

shopt -s expand_aliases
alias fastp="apptainer exec /usr/local/biotools/f/fastp:0.23.4--h5f740d0_0 fastp"

[ ! -e ${clean_data} ] && mkdir ${clean_data}
[ ! -e out/qc ] && mkdir -p out/qc

paired_fastaq=($(ls ${raw_data}/${sample} | sort -V))
paired_1=${paired_fastaq[0]}
paired_2=${paired_fastaq[1]}

fastp \
  -i ${raw_data}/${sample}/${paired_1} \
  -I ${raw_data}/${sample}/${paired_2} \
  -o ${clean_data}/${sample}_clean_1.fq.gz \
  -O ${clean_data}/${sample}_clean_2.fq.gz \
  -h out/qc/${sample}_qc.html \
  -q 30 -u 30

rm -r ${raw_data}/${sample}
