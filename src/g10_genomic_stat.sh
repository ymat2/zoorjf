#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"

proj=~/RJF
vcf=${proj}/vcf/RJF.snp.vcf.gz
workdir=${proj}/gscan

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}
mkdir pi tajimasD pop

## Per individual heterozygosity
# vcftools --gzvcf ${vcf} --het --maf 0.01 --max-missing 0.1 --out het/RJF

## make population text
for pop in $(cat ../treemix/treemix.list); do
  cat ../treemix/treemix.clust | grep ${pop} | awk '{ print $1 }' > pop/${pop}.txt
done

## Pi and Tajima's D
for pop in $(cat ../treemix/treemix.list); do
  vcftools --gzvcf ${vcf} --mac 2 --window-pi 50000 --window-pi-step 50000 --keep pop/${pop}.txt --out pi/${pop}
  vcftools --gzvcf ${vcf} --mac 2 --TajimaD 50000 --keep pop/${pop}.txt --out tajimasD/${pop}
done
