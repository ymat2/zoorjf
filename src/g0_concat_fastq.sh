#$ -S /bin/bash
#$ -cwd
#$ -l s_vmem=16G
#$ -l mem_req=16G
#$ -o /dev/null
#$ -e /dev/null


# raw_data1=~/raw_data/RJF_ABRC
raw_data2=~/raw_data/RJF_KMT
raw_data3=~/raw_data/RJF_iZoo
proj=~/RJF

cd ${proj}

for d in ${raw_data2} ${raw_data3}; do
  samples=($(ls ${d}))
  for sample in ${samples[@]}; do
    [ ! -e uncleaned/${sample} ] && mkdir -p uncleaned/${sample}
    paired_1=($(ls ${d}/${sample}/*_1.fq.gz))
    paired_2=($(ls ${d}/${sample}/*_2.fq.gz))
    cat ${paired_1[@]} > uncleaned/${sample}/${sample}_1.fq.gz
    cat ${paired_2[@]} > uncleaned/${sample}/${sample}_2.fq.gz
  done
done
