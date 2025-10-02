#!/usr/bin/env bash
# use trim5p and/or trim3p to remove artifactual patterns of mismatches at the ends of the reads
# examples trim 15bp for 5p only
# -double will only use doubled sequenced part for estimation
# double seems to make reads more likely be new
# make sure having no4sU sample
gedi -e Slam \
-progress -D -trim5p 15 -nthreads 22 \
-snpConv 0.01 -snppval 0.1 \
-modelall \
-no4sUpattern BMP_seth \
-genomic gedi_homo \
-prefix gedi_out/BMP_mTe \
-reads cit.bamlist.cit
