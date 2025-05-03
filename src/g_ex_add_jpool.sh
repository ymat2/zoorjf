#!/bin/bash

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias plink1="apptainer exec /usr/local/biotools/p/plink:1.90b6.21--hec16e2b_4 plink"

proj=~/RJF
workdir=${proj}/_add_jpool
vcf=RJF_jadd.snp.vcf
prefix=RJF_jadd.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

bcftools merge ../vcf/RJF.snp.vcf.gz ~/RJF2/vcf/RJF2.snp.vcf.gz -0 -Oz > ${vcf}
bcftools index ${vcf}

bcftools filter ${proj}/vcf/RJF.snp.vcf.gz -i "DP<125" -Oz -o ${prefix}.filter.vcf.gz
bcftools index ${prefix}.filter.vcf.gz

plink1 --vcf ${proj}/vcf/RJF.snp.vcf.gz \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --maf 0.01 \
  --geno 0.1 \
  --indep-pairwise 150 50 0.2 \
  --out ${prefix}

plink1 --vcf ${proj}/vcf/RJF.snp.vcf.gz \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --extract ${prefix}.prune.in \
  --make-bed \
  --out ${prefix}

rm ${prefix}.prune*

plink1 --bfile ${prefix} \
  --pca header tabs \
  --allow-extra-chr \
  --double-id \
  --out ${prefix}.pca
