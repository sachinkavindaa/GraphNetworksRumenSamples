#!/bin/bash

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kofamscan_env

# Set paths
KOFAMSCAN_DIR="/work/samodha/sachin/kofam_scan"
INPUT_DIR="/work/samodha/sachin/GraphN/split_fasta"
OUTPUT_DIR="/work/samodha/sachin/GraphN/kofam_results"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all FASTA files in the input directory
for fasta in "$INPUT_DIR"/*.faa; do
    filename=$(basename "$fasta" .faa)
    echo "Processing $filename..."

    "$KOFAMSCAN_DIR"/exec_annotation \
        -o "$OUTPUT_DIR/${filename}_ko.tsv" \
        -p "$KOFAMSCAN_DIR/profiles" \
        -k "$KOFAMSCAN_DIR/ko_list" \
        --cpu 4 \
        "$fasta"
done

echo "KO assignment completed!"
