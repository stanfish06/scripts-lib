#!/usr/bin/env bash

# use this in case you got many polyA and polyT contamination (e.g. >1% in RNAseq)
module load Bioinformatics
module load cutadapt
module load fastqc

echo begin
cutadapt \
	-j 12 \
	-a A{20} -a T{20} \
	-A A{20} -A T{20} \
	-m 50 \
	-o BC3_1.trim_poly.fq.gz -p BC3_2.trim_poly.fq.gz \
	./BC3_trimmed_1_val_1.fq.gz ./BC3_trimmed_2_val_2.fq.gz
echo start qc
fastqc -t 2 BC3_1.trim_poly.fq.gz BC3_2.trim_poly.fq.gz
echo done
