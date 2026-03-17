#!/usr/bin/env bash
module load use.own
module load gedi
RUN_CIT_PREP=true
RUN_GENOME_PREP=false
RUN_SLAM=false

# Generate CIT file
if $RUN_CIT_PREP; then
echo "Step 1/3: Generate CIT file for the bam list"
cat > crosstalk.bamlist <<EOF
star_salmon/A_1hr.sorted.bam
star_salmon/A_8hr.sorted.bam
star_salmon/A_LDN_1hr_2.sorted.bam
star_salmon/A_LDN_8hr.sorted.bam
star_salmon/B_1hr.sorted.bam
star_salmon/B_8hr.sorted.bam
star_salmon/B_A_1hr.sorted.bam
star_salmon/B_A_8hr.sorted.bam
star_salmon/B_A_SB_1hr_2.sorted.bam
star_salmon/B_A_SB_8hr.sorted.bam
star_salmon/B_SB_1hr.sorted.bam
star_salmon/B_SB_8hr.sorted.bam
star_salmon/ctrl_lb1_8hr.sorted.bam
star_salmon/ctrl_lb_1hr.sorted.bam
star_salmon/ctrl_lb2_8hr.sorted.bam
star_salmon/ctrl_unl_1hr.sorted.bam
star_salmon/ctrl_unl_8hr.sorted.bam
star_salmon/GDF11_1hr.sorted.bam
star_salmon/GDF11_8hr.sorted.bam
star_salmon/SB_1hr.sorted.bam
star_salmon/SB_8hr.sorted.bam
EOF
bamlist2cit crosstalk.bamlist
fi

# Prep genome
if $RUN_GENOME_PREP; then
SOURCE="/nfs/turbo/umms-iheemske/reference-genomes/human/GRCh38/GENCODE/v49"
FASTA="GRCh38.p14.primary_assembly.genome.fa"
GTF="gencode.v49.primary_assembly.annotation.gtf"
echo "Step 2/3: generate genome index"
cp "${SOURCE}/${FASTA}.gz" . && gunzip -f "${FASTA}.gz"
cp "${SOURCE}/${GTF}.gz" . && gunzip -f "${GTF}.gz"
gedi -e IndexGenome \
-s $FASTA \
-a $GTF \
-n homo_hg38_gedi
fi

if $RUN_SLAM; then
# use trim5p and/or trim3p to remove artifactual patterns of mismatches at the ends of the reads
# examples trim 15bp for 5p only
# -double will only use doubled sequenced part for estimation
# double seems to make reads more likely be new
# make sure having no4sU sample
echo "Step 3/3: Run SLAM"
gedi -e Slam \
-progress -D -trim5p 15 -nthreads 22 \
-snpConv 0.01 -snppval 0.1 \
-modelall \
-no4sUpattern unl \
-genomic homo_hg38_gedi \
-prefix gedi_out/BMP_ACT_crosstalk \
-reads cit.bamlist.cit
fi
