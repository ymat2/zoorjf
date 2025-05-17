#!/bin/bash

#SBATCH -a 1-130
#SBATCH --mem 32G


samples=($(ls bam | sort -V))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

workdir=~/RJF/psmc
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

#psmc psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${sample}/${sample}.4.psmc ${sample}/${sample}.split.filtered.psmcfa
#psmc psmc -N25 -t15 -r5 -p "2+2+25*2+4+6" -o ${sample}/${sample}.22.psmc ${sample}/${sample}.split.filtered.psmcfa
psmc psmc -N25 -t15 -r5 -p "1+1+1+1+25*2+4+6" -o ${sample}/${sample}.1111.psmc ${sample}/${sample}.split.filtered.psmcfa
rm ${sample}/${sample}.split.filtered.psmcfa
