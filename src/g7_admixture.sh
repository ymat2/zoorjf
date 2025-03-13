#!/bin/bash

#SBATCH -a 1-10
#SBATCH --mem 64G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

seq_ids=(1 2 3 4 5 6 7 8 9 10)
seq_id=${seq_ids[$SLURM_ARRAY_TASK_ID-1]}

shopt -s expand_aliases
alias admixture="apptainer exec /usr/local/biotools/a/admixture:1.3.0--0 admixture"

workdir=~/RJF/structure
prefix=RJF.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

echo Starting ADMIXTURE K: ${seq_id}
admixture --cv ${prefix}.bed ${seq_id} | tee log${seq_id}.out

grep -h CV log*.out > CV-error.txt
