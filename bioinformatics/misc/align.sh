#!/usr/bin/env bash

module load Bioinformatics
module load star/2.7.11a-hdp2onj

STAR \
	--genomeDir ./star \
	--readFilesIn BC3_1.trim_poly.fq.gz BC3_2.trim_poly.fq.gz \
	--readFilesCommand zcat \
	--runThreadN 12 \
	--alignEndsType Local \
	--outFileNamePrefix BC3_realign_ \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts
