#!/usr/bin/bash
awk 'BEGIN{FS=OFS="\t"} !/^#/{
  gsub(/=/," "); gsub(/;/,";");
  split($9,a,";"); $9="";
  for(i in a){ if(a[i]!=""){ split(a[i],b," "); $9=$9 sprintf("%s \"%s\"; ", b[1], b[2]) } }
  print
}' ./Homo_sapiens.GRCh38.regulatory_features.v115.gff3 > ./HS_reg.gtf

