#!/bin/bash

#SBATCH -a 1-17
#SBATCH --time 60
#SBATCH --partition short

shopt -s expand_aliases
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"

proj=~/RJF
vcf=${proj}/vcf/RJF.snp.vcf.gz
workdir=${proj}/gscan

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

pops=($(cat ../treemix/treemix.list))
pop=${pops[$SLURM_ARRAY_TASK_ID-1]}

## make population text
echo ${pop}
cat ../treemix/treemix.clust | grep ${pop} | awk '{ print $1 }' > pop/${pop}.txt

## Pi and Tajima's D
vcftools --gzvcf ${vcf} --window-pi 50000 --keep pop/${pop}.txt --out pi/${pop}
vcftools --gzvcf ${vcf} --TajimaD 50000 --keep pop/${pop}.txt --out tajimasD/${pop}

echo Done!!
