#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias plink1="apptainer exec /usr/local/biotools/p/plink:1.90b6.21--hec16e2b_4 plink"

workdir=~/RJF/structure
prefix=RJF.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

plink1 --bfile ${prefix} --pca --allow-extra-chr --double-id --out ${prefix}.pca

# For ADMIXTURE
[ ! -e ~/RJF/admixture ] && mkdir ~/RJF/admixture
cp ${prefix}.bim ~/RJF/admixture/
cp ${prefix}.bed ~/RJF/admixture/
cp ${prefix}.fam ~/RJF/admixture/
sed -i -e 's/NW_//g' ~/RJF/admixture/${prefix}.bim
sed -i -e 's/NC_//g' ~/RJF/admixture/${prefix}.bim
sed -i -e 's/\.1//g' ~/RJF/admixture/${prefix}.bim
