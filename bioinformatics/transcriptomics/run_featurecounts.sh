#!/bin/bash
# Run featureCounts with gene_id grouping for paired-end data
module load Bioinformatics subread

WORKDIR="/gpfs/accounts/iheemske_root/iheemske0/zyyu/260120-Seth/star_salmon"
GTF="${WORKDIR}/gencode.v49.primary_assembly.annotation.filtered.gtf"
OUTPUT="${WORKDIR}/featurecounts.merged.gene_counts.tsv"

cd "$WORKDIR"

BAMS=(
  A_1hr.markdup.sorted.bam
  A_8hr.markdup.sorted.bam
  A_LDN_1hr_2.markdup.sorted.bam
  A_LDN_8hr.markdup.sorted.bam
  B_1hr.markdup.sorted.bam
  B_8hr.markdup.sorted.bam
  B_A_1hr.markdup.sorted.bam
  B_A_8hr.markdup.sorted.bam
  B_A_SB_1hr_2.markdup.sorted.bam
  B_A_SB_8hr.markdup.sorted.bam
  B_SB_1hr.markdup.sorted.bam
  B_SB_8hr.markdup.sorted.bam
  ctrl_lb_1hr.markdup.sorted.bam
  ctrl_lb1_8hr.markdup.sorted.bam
  ctrl_lb2_8hr.markdup.sorted.bam
  ctrl_unl_1hr.markdup.sorted.bam
  ctrl_unl_8hr.markdup.sorted.bam
  GDF11_1hr.markdup.sorted.bam
  GDF11_8hr.markdup.sorted.bam
  SB_1hr.markdup.sorted.bam
  SB_8hr.markdup.sorted.bam
)

# Run featureCounts with paired-end mode
featureCounts \
  -p --countReadPairs \
  -B -C \
  -g gene_id \
  -t exon \
  -T 24 \
  -s 0 \
  -a "$GTF" \
  -o "$OUTPUT" \
  "${BAMS[@]}"

# Clean up output: remove comment line, extract gene_id + counts, simplify column names
tail -n +2 "$OUTPUT" | \
  cut -f1,7- | \
  sed '1s/\.markdup\.sorted\.bam//g' > "${WORKDIR}/featurecounts.merged.gene_counts.clean.tsv"

echo "output: featurecounts.merged.gene_counts.clean.tsv"
