#!/bin/bash

shopt -s expand_aliases
alias plink1="apptainer exec /usr/local/biotools/p/plink:1.90b6.21--hec16e2b_4 plink"

workdir=~/RJF/structure
prefix=RJF.snp

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

plink1 --bfile ${prefix} \
  --pca header tabs var-wts \
  --allow-extra-chr \
  --double-id \
  --out ${prefix}.pca

cat ${prefix}.pca.eigenvec.var | awk '{print $2 "\t" $5}' > ${prefix}.pca.PC1.var
rm ${prefix}.pca.eigenvec.var

# For ADMIXTURE
[ ! -e ~/RJF/admixture ] && mkdir ~/RJF/admixture
cp ${prefix}.bim ~/RJF/admixture/
cp ${prefix}.bed ~/RJF/admixture/
cp ${prefix}.fam ~/RJF/admixture/
sed -i -e 's/NW_//g' ~/RJF/admixture/${prefix}.bim
sed -i -e 's/NC_//g' ~/RJF/admixture/${prefix}.bim
sed -i -e 's/\.1//g' ~/RJF/admixture/${prefix}.bim
