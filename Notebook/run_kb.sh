#!/bin/bash
FASTQ_DIR="/home/nakagawa/datasets/SRR_cortex/SRR_cortex_scRNA"
OUT_DIR="/home/nakagawa/datasets/SRR_cortex/SRR_cortex_scRNA/processed_h5ad"
INDEX="/home/nakagawa/datasets/genome/mm10/transcriptome.idx"
T2G="/home/nakagawa/datasets/genome/mm10/transcripts_to_genes.txt"

mkdir -p "$OUT_DIR"

for R2 in ${FASTQ_DIR}/*_2.fastq.gz; do
    SAMPLE=$(basename "$R2" _2.fastq.gz)
    echo "---------------------------------------------------"
    echo "Processing $SAMPLE..."
    SAMPLE_OUT="$OUT_DIR/$SAMPLE"
    R1="$FASTQ_DIR/${SAMPLE}_1.fastq.gz"
    if [[ ! -f "$R1" ]]; then echo "Missing R1 for $SAMPLE"; continue; fi
    rm -rf "$SAMPLE_OUT"
    echo "Trying 10xv3..."
    kb count -i "$INDEX" -g "$T2G" -x 10xv3 -o "$SAMPLE_OUT" --h5ad "$R1" "$R2"
    STATUS=$?
    if [[ $STATUS -ne 0 ]]; then
        echo "Retrying with 10xv2..."
        rm -rf "$SAMPLE_OUT"
        kb count -i "$INDEX" -g "$T2G" -x 10xv2 -o "$SAMPLE_OUT" --h5ad "$R1" "$R2"
        STATUS=$?
    fi
    if [[ $STATUS -eq 0 ]]; then echo "Finished $SAMPLE successfully"; else echo "FAILED: $SAMPLE"; fi
done
