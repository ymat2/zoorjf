#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias plink1="apptainer exec /usr/local/biotools/p/plink:1.90b6.21--hec16e2b_4 plink"
alias plink2="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"

proj=~/RJF
workdir=${proj}/roh
vcf=${proj}/vcf/RJF.snp.vcf.gz
prefix=RJF.snp
maf=0.01
geno=0.1

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

plink1 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --geno ${gneo} \
  --maf ${maf} \
  --not-chr NC_052571.1, NC_052572.1, NC_053523.1 \
  --set-missing-var-ids @:# \
  --homozyg-snp 50 \
  --homozyg-window-snp 50 \
  --homozyg-window-missing 5 \
  --homozyg-window-het 3 \
  --homozyg-kb 300 \
  --homozyg-density 50 \
  --homozyg-gap 1000 \
  --homozyg-window-threshold 0.05 \
  --out ${prefix}.roh


plink1 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --geno 0 \
  --maf ${maf} \
  --het \
  --not-chr NC_052571.1, NC_052572.1, NC_053523.1 \
  --set-missing-var-ids @:# \
  --out ${prefix}
