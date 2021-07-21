#!/bin/bash
# This is just a crude way of getting sample names from single-sample GVCFs

cat gvcfs_to_merge.txt | while read gvcf
do
  extension="${gvcf##*.}"

  if [[ "$extension" == "gvcf" || "$extension" == "vcf" ]]
  then
    gsutil cat ${gvcf} | grep -m1 "#CHROM" | awk '{print $10}' >> sample_names.txt

  elif [[ "$extension" == "gz" ]]
  then
    gsutil cat ${gvcf} | gunzip | grep -m1 "#CHROM" | awk '{print $10}' >> sample_names.txt

  else
    usage
    stop "unrecognised extension $extension"
  fi
done