#$ -S /bin/bash
#$ -cwd
#$ -t 1-16:1
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -o /dev/null
#$ -e /dev/null


samples=($(awk -F '\t' 'NR>1 {print $1}' data/sra_accession_add.tsv))
sample=${samples[$SGE_TASK_ID-1]}
reference=~/ref/grcg7b/GCF_016699485.2.fa

shopt -s expand_aliases
alias bwa="apptainer exec /usr/local/biotools/b/bwa:0.7.17--h5bf99c6_8 bwa"
alias samtools="apptainer exec /usr/local/biotools/s/samtools:1.18--h50ea8bc_1 samtools"

mkdir -p bam/${sample}

bwa mem -t 4 -M ${reference} clean/${sample}_clean_1.fq.gz clean/${sample}_clean_2.fq.gz |
  samtools view -b -F 4 > bam/${sample}/${sample}.bam
samtools collate -Ou bam/${sample}/${sample}.bam | \
  samtools fixmate - - -mu | \
  samtools sort - -u | \
  samtools markdup - bam/${sample}/${sample}.cfsm.bam
samtools index bam/${sample}/${sample}.cfsm.bam

samtools flagstat -O tsv bam/${sample}/${sample}.bam > bam/${sample}/${sample}.stat
samtools flagstat -O tsv bam/${sample}/${sample}.cfsm.bam > bam/${sample}/${sample}.cfsm.stat
samtools coverage bam/${sample}/${sample}.cfsm.bam > bam/${sample}/${sample}.cov

rm bam/${sample}/${sample}.bam
