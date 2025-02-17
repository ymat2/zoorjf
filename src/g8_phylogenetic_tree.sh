#$ -S /bin/bash
#$ -cwd
#$ -l medium
#$ -l mem_req=512G
#$ -l s_vmem=512G
#$ -o /dev/null
#$ -e /dev/null


shopt -s expand_aliases
alias plink2="apptainer exec /usr/local/biotools/p/plink2:2.00a5--h4ac6f70_0 plink2"
alias rapidnj="apptainer exec /usr/local/biotools/r/rapidnj:2.3.2--h4ac6f70_4 rapidnj"

workdir=~/RJF/phylo/
vcf=~/RJF/vcf/RJF.snp.vcf.gz
prefix=RJF.sub

[ ! -e ${workdir} ] && mkdir -p ${workdir}
cd ${workdir}

plink2 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --maf 0.01 \
  --not-chr NC_052571.1, NC_052572.1, NC_053523.1 \
  --indep-pairwise 50 10 0.2 \
  --out ${prefix}

plink2 --vcf ${vcf} \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids @:# \
  --extract ${prefix}.prune.in \
  --export vcf id-paste=iid \
  --out ${prefix}

/usr/bin/python3 ~/bin/vcf2phylip/vcf2phylip.py -i ${prefix}.vcf
rm ${prefix}.vcf* ${prefix}.prune*

poppy alnkit --mode trim -i ${prefix}.min4.phy -o ${prefix}.varsites.fa --format fasta
rm ${prefix}.min4.phy

rapidnj ${prefix}.varsites.fa -i fa -o t -b 100 --evolution-model kim -x ${prefix}.nj.tree
