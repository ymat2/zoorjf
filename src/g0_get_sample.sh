#!/bin/bash

#SBATCH -a 1-16
#SBATCH --mem 32G
#SBATCH -o /dev/null
#SBATCH -e /dev/null


sras=($(awk -F '\t' 'NR>1 {print $1}' data/sra_accession.tsv))
sra=${sras[$SLURM_ARRAY_TASK_ID-1]}

workdir=~/RJF/sra
uncleaned=~/RJF/uncleaned/

[ ! -e ${workdir} ] && mkdir -p ${workdir}
[ ! -e ${uncleaned} ] && mkdir -p ${uncleaned}

cd ${workdir}

prefetch --max-size 50G ${sra}
fasterq-dump ${sra}/${sra}.sra --outdir ${uncleaned}/${sra}
rm -r ${sra}
