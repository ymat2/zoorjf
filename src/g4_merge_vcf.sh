#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

proj=~/RJF
cd ${proj}

poppy summary -i ~/RJF/bam -o ~/RJF/out/flag_summary.tsv --mode flag --suffix stat
poppy summary -i ~/RJF/bam -o ~/RJF/out/mapping_summary.tsv --mode bam --suffix cov
#poppy summary -i ~/RJF/bam -o ~/RJF/out/vcf_summary.tsv --mode vcf --suffix vcf.stat

samples=($(ls ~/RJF/bam | sort -V))
for s in ${samples[@]}; do echo ~/RJF/bam/$s/$s.q.vcf.gz >> data/use_samples.txt ; done
cat ~/RJF/out/mapping_summary.tsv | awk '{ if ($2 < 3 || $3 < 90) { print $1 }}' | sort
# edit data/use_samples.txt

workdir=./vcf
vcf=${workdir}/RJF.vcf.gz
snponly=${workdir}/RJF.snp.vcf.gz

[ ! -e ${workdir} ] && mkdir ${workdir}

samples=($(cat data/use_samples.txt | grep -v -e "^#"))
bcftools merge ${samples[@]} -0 -Oz > ${vcf}
bcftools index ${vcf}

bcftools view -v snps -m2 -M2 -Oz ${vcf} > ${snponly}
bcftools index ${snponly}

## Basic statistics
bcftools stats --depth 500,1000,1 -s - ${vcf} > out/RJF.vcf.stats
