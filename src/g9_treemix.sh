#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias plink1="apptainer exec /usr/local/biotools/p/plink:1.90b6.21--hec16e2b_4 plink"
alias treemix="apptainer exec /usr/local/biotools/t/treemix:1.13--h125836d_9 treemix"

workdir=~/RJF/treemix
vcf=~/RJF/vcf/RJF.snp.vcf.gz
bfile=~/RJF/structure/RJF.snp
prefix=RJF.sub

[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

## Make .clust file: make_treemix_clust.R
awk '{print $3}' treemix.clust | sort | uniq > treemix.list

## Make Treemix input file

plink1 --bfile ${bfile} \
  --freq \
  --missing \
  --within treemix.clust \
  --not-chr NC_052571.1, NC_052572.1, NC_053523.1 \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --out ${prefix}

# gzip ${prefix}.frq.strat

python3 ../src/plink2treemix.py -i ${prefix}.frq.strat -o ${prefix}.treemix.frq
rm ${prefix}.*miss ${prefix}.frq.strat ${prefix}.nosex
gzip ${prefix}.treemix.frq

for m in {1..8}; do
  for i in {1..5}; do
    echo Running Treemix: M=${m}, iteration=${i} ...
    treemix -i ${prefix}.treemix.frq.gz -m ${m} -o treemix.${i}.${m} -bootstrap -k 500 -root gallus -seed ${i} > treemix_${i}.${m}_log
  done
done
