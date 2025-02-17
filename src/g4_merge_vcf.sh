#$ -S /bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null


shopt -s expand_aliases
alias bcftools="apptainer exec /usr/local/biotools/b/bcftools:1.18--h8b25389_0 bcftools"

proj=~/RJF
cd ${proj}

poppy summary -i ~/RJF/bam -o ~/RJF/out/flag_summary.tsv --mode flag --suffix stat
poppy summary -i ~/RJF/bam -o ~/RJF/out/mapping_summary.tsv --mode bam --suffix cov
poppy summary -i ~/RJF/bam -o ~/RJF/out/vcf_summary.tsv --mode vcf --suffix vcf.stat

echo -e "sample\tn_snp\tn_hom\tn_het\tn_mnp\tn_indel\tn_missing\tn_multi" > out/snp_count.txt
for s in $(ls bam); do
  line=$(poppy vcfkit --mode count -i ~/RJF/bam/${s}/${s}.q.vcf.gz | grep "bam")
  echo ${line} >> out/snp_count.txt
done

samples=($(ls ~/RJF/bam | sort -V))
for s in ${samples[@]}; do echo ~/RJF/bam/$s/$s.q.vcf.gz >> data/use_samples.txt ; done
# edit data/use_samples.txt

workdir=./vcf
vcf=${workdir}/RJF.vcf.gz

[ ! -e ${workdir} ] && mkdir ${workdir}

samples=($(cat data/use_samples.txt | grep -v -e "^#"))
bcftools merge ${samples[@]} -0 -Oz > ${vcf}
bcftools index ${vcf}
