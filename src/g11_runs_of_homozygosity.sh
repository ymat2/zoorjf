#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias plink1="apptainer exec /usr/local/biotools/p/plink:1.90b6.21--hec16e2b_4 plink"

proj=~/RJF
workdir=${proj}/roh
snponly=${proj}/vcf/RJF.snp.vcf.gz
prefix=RJF.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

plink1 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --maf 0.01 \
  --not-chr NC_053523.1 \
  --indep-pairwise 150 50 0.2 \
  --out ${prefix}

plink1 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --extract ${prefix}.prune.in \
  --make-bed \
  --out ${prefix}

rm ${prefix}.prune.in ${prefix}.prune.out

plink1 --bfile ${prefix} \
  --allow-extra-chr \
  --double-id \
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
