#!/bin/bash

#SBATCH -a 1-10
#SBATCH --mem 32G


sample=TMP0001M

workdir=~/RJF/psmc_param
reference=~/ref/grcg7b/GCF_016699485.2.fa

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias psmc="apptainer exec /usr/local/biotools/p/psmc:0.6.5--h43eeafb_2"

[ ! -e ${workdir}/${sample} ] && mkdir -p ${workdir}/${sample}
cd ${workdir}

bcftools consensus -f ${reference} ~/RJF/bam/${sample}/${sample}.q.vcf.gz -I | gzip > ${sample}/${sample}.consensus.fa.gz
psmc fq2psmcfa -s 50 ${sample}/${sample}.consensus.fa.gz > ${sample}/${sample}.consensus.psmcfa
psmc splitfa ${sample}/${sample}.consensus.psmcfa > ${sample}/${sample}.split.filtered.psmcfa
rm ${sample}/${sample}.consensus.fa.gz ${sample}/${sample}.consensus.psmcfa

psmc psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${sample}/${sample}.4.psmc ${sample}/${sample}.split.filtered.psmcfa
psmc psmc -N25 -t15 -r5 -p "2+2+25*2+4+6" -o ${sample}/${sample}.22.psmc ${sample}/${sample}.split.filtered.psmcfa
psmc psmc -N25 -t15 -r5 -p "1+1+1+1+25*2+4+6" -o ${sample}/${sample}.1111.psmc ${sample}/${sample}.split.filtered.psmcfa
rm ${sample}/${sample}.split.filtered.psmcfa

#perl ../src/psmc_plot.pl -s 50 -X 50000000 -p -g 1 -R -x 1000 -u 1.91e-9 ${sample}/${sample}.plot ${sample}/${sample}.22.psmc
#perl ../src/psmc_plot.pl -s 50 -X 50000000 -p -g 1 -R -x 1000 -u 1.91e-9 ${sample}/${sample}.plot ${sample}/${sample}.1111.psmc
