#!/bin/bash

#SBATCH -o /dev/null
#SBATCH -e /dev/null

shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"
alias bedtools="apptainer exec /usr/local/biotools/b/bedtools:2.31.0--h468198e_0 bedtools"
alias rapidnj="apptainer exec /usr/local/biotools/r/rapidnj:2.3.2--h4ac6f70_4 rapidnj"

proj=~/RJF
workdir=${proj}/fst
vcf=${proj}/vcf/RJF.snp.vcf.gz

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}

bcftools query -l ${vcf} | grep -E "HoU|RJF|TMP" > zoo.txt
bcftools query -l ${vcf} | grep -E -v "HoU|RJF|TMP" > other.txt

vcftools --gzvcf ${vcf} --weir-fst-pop zoo.txt --weir-fst-pop other.txt \
  --fst-window-size 50000 --fst-window-step 50000 --out zoo_vs_other


### phylogenetic analysis of high MEAN_FST region

bcftools view --regions NC_052533.1:52550001-53050000 ${vcf} > wild_vs_zoo.highfst.vcf
/usr/bin/python3 ~/bin/vcf2phylip/vcf2phylip.py -i wild_vs_zoo.highfst.vcf
rm wild_vs_zoo.highfst.vcf*

poppy alnkit --mode trim -i wild_vs_zoo.highfst.min4.phy -o wild_vs_zoo.highfst.varsites.fa --format fasta
rapidnj wild_vs_zoo.highfst.varsites.fa -i fa -o t -b 100 --evolution-model kim -x wild_vs_zoo.highfst.nj.tree
rm wild_vs_zoo.highfst.min4.phy wild_vs_zoo.highfst.varsites.fa
