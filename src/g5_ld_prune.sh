#$ -S /bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias plink1="apptainer exec /usr/local/biotools/p/plink:1.90b6.21--hec16e2b_4 plink"

proj=~/RJF
workdir=${proj}/structure
vcf=${proj}/vcf/RJF.vcf.gz
snponly=${proj}/vcf/RJF.snp.vcf.gz
prefix=RJF.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

bcftools view -v snps -m2 -M2 -Oz ${vcf} > ${snponly}
bcftools index ${snponly}

plink1 --vcf ${snponly} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --maf 0.01 \
  --indep-pairwise 50 10 0.2 \
  --out ${prefix}

plink1 --vcf ${snponly} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --extract ${prefix}.prune.in \
  --make-bed \
  --out ${prefix}

rm ${prefix}.prune*
