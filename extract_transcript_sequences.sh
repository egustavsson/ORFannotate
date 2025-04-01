#!/bin/bash

# Check for required arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input.sam> <annotations.gtf> <output.bam> <output.fasta>"
    exit 1
fi

# Assign input arguments
SAM_FILE=$1
GTF_FILE=$2
OUTPUT_BAM=$3
OUTPUT_FASTA=$4

# Generate temporary filenames
SORTED_BAM="sorted.bam"
BED_FILE="transcripts.bed"

echo "Converting SAM to sorted BAM..."
samtools view -bS "$SAM_FILE" | samtools sort -o "$SORTED_BAM"
samtools index "$SORTED_BAM"

echo "Extracting transcript coordinates from GTF..."
awk '$3 == "transcript" {print $1, $4, $5, $9}' "$GTF_FILE" > "$BED_FILE"

echo "Filtering BAM file to retain transcript-aligned reads..."
samtools view -L "$BED_FILE" -b "$SORTED_BAM" > "$OUTPUT_BAM"

echo "Converting BAM to FASTA..."
samtools fasta "$OUTPUT_BAM" > "$OUTPUT_FASTA"

echo "Transcript sequences saved to $OUTPUT_FASTA."
