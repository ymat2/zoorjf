#!/bin/bash

samples=($(ls bam | grep -E "HoU|IZ|RJF|TMP"))
for s in ${samples[@]}; do
  cat bam/${s}/${s}.cov | grep -E "NC_052571.1|NC_052572.1" | awk -v s=${s} '{print s"\t"$1"\t"$4"\t"$6"\t"$7}'
done > out/sex_chrom_cov.tsv
