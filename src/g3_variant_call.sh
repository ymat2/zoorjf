#$ -S /bin/bash
#$ -cwd
#$ -t 1-16
#$ -l s_vmem=63G
#$ -l mem_req=63G
#$ -o /dev/null
#$ -e /dev/null

samples=($(ls bam | sort -V))
sample=${samples[$SGE_TASK_ID-1]}
reference=~/ref/grcg7b/GCF_016699485.2.fa

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

ave_coverage=$(awk '{NR>1 && total_dp += $5*$7; total_cb += $5} END {print total_dp/total_cb}' bam/${sample}/${sample}.cov)

bcftools mpileup -f ${reference} bam/${sample}/${sample}.cfsm.bam |
  bcftools call -vm |
  bcftools filter -i "QUAL>30 && DP>${ave_coverage}/2 && DP<${ave_coverage}*2" --set-GTs . -Oz -o bam/${sample}/${sample}.q.vcf.gz
bcftools index bam/${sample}/${sample}.q.vcf.gz

bcftools stats bam/${sample}/${sample}.q.vcf.gz > bam/${sample}/${sample}.vcf.stat
