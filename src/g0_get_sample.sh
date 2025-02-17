#$ -S /bin/bash
#$ -cwd
#$ -t 1-16:1
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -o /dev/null
#$ -e /dev/null


sras=($(awk -F '\t' 'NR>1 {print $1}' data/sra_accession.tsv))
sra=${sras[$SGE_TASK_ID-1]}

workdir=~/RJF/sra
uncleaned=~/RJF/uncleaned/

[ ! -e ${workdir} ] && mkdir -p ${workdir}
[ ! -e ${uncleaned} ] && mkdir -p ${uncleaned}

cd ${workdir}

prefetch --max-size 50G ${sra}
fasterq-dump ${sra}/${sra}.sra --outdir ${uncleaned}/${sra}
rm -r ${sra}
