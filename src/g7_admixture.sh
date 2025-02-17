#$ -S /bin/bash
#$ -cwd
#$ -t 1-8:1
#$ -l mem_req=63G
#$ -l s_vmem=63G
#$ -o /dev/null
#$ -e /dev/null

seq_ids=(1 2 3 4 5 6 7 8)
seq_id=${seq_ids[$SGE_TASK_ID-1]}

shopt -s expand_aliases
alias admixture="apptainer exec /usr/local/biotools/a/admixture:1.3.0--0 admixture"

workdir=~/RJF/structure
prefix=RJF.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

echo Starting ADMIXTURE K: ${seq_id}
admixture --cv ${prefix}.bed ${seq_id} | tee log${seq_id}.out

grep -h CV log*.out > CV-error.txt
