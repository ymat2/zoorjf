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
vcftools --gzvcf ${vcf} --het --maf 0.01 --max-missing 0.1 --out het/RJF

## Windonwed Pi and Tajima's D across populations
vcftools --gzvcf ${vcf} --window-pi 50000 --out pi/RJF_all
vcftools --gzvcf ${vcf} --TajimaD 50000 --out tajimasD/RJF_all
