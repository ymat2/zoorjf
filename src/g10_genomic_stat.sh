#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -o /dev/null
#$ -e /dev/null


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"
alias vcftools="apptainer exec /usr/local/biotools/v/vcftools:0.1.16--h9a82719_5 vcftools"

proj=~/RJF
vcf=${proj}/vcf/RJF.snp.vcf.gz
workdir=${proj}/gscan

[ ! -e ${workdir} ] && mkdir ${workdir}
cd ${workdir}
mkdir pi tajimasD het

## Per individual heterozygosity
vcftools --gzvcf ${vcf} --het --out het/RJF

## Pi
for pop in GGg GGj GGsc GGss GGst guax idn izoo kmt thai viet yunn LDH BLH; do
  vcftools --gzvcf ${vcf} --mac 2 --window-pi 10000 --window-pi-step 10000 --keep ${pop}.txt --out pi/${pop}
done

## Tajima's D
for pop in GGg GGj GGsc GGss GGst guax idn izoo kmt thai viet yunn LDH BLH; do
  vcftools --gzvcf ${vcf} --mac 2 --TajimaD 10000 --keep ${pop}.txt --out tajimasD/${pop}
done
