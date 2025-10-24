#!/bin/bash

REG_GTF="./HS_reg.gtf"
THREADS=24
BAM_FILES=()
for sample in SRR24163129 SRR24163130 SRR24163131 SRR24163132 SRR24163133 SRR24163134; do
    BAM_FILE="${sample}Aligned.sortedByCoord.out.bam"
    if [ -f "$BAM_FILE" ]; then
        BAM_FILES+=("$BAM_FILE")
    else
        echo "Warning: $BAM_FILE not found, skipping."
    fi
done

featureCounts \
    -F GFF \
    -t enhancer \
    -g ID \
    -T $THREADS \
    -p \
    -O \
    --fracOverlap 0.2 \
    -a "$REG_GTF" \
    -o all_samples_enhancer_counts.txt \
    "${BAM_FILES[@]}"

