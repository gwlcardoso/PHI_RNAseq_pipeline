#!/usr/bin/env bash
# ========================
# Gene counting with featureCounts
# Author: Guilherme Wenceslau
# Year: 2025
# ========================

set -euo pipefail
IFS=$'\n\t'
shopt -s nullglob

# ------------------------
# User configuration
# ------------------------
THREADS=8
STRANDNESS=0
FEATURE_TYPE="CDS"
ATTRIBUTE_TYPE="transcript_id"

PROJECT_DIR="/home/gwlcardoso/Projects/phi_rnaseq"

# Input annotation files
HOST_GFF="${PROJECT_DIR}/nextseq_data/ref_gen/PA14.gff"
PHAGE_GFF="${PROJECT_DIR}/nextseq_data/ref_gen/ZC01.gff"

# Output annotation files
ANNOTATION_DIR="${PROJECT_DIR}/nextseq_data/ref_gen/combined_annotation"

COMBINED_GFF="${ANNOTATION_DIR}/PA14_ZC01.gff"
COMBINED_GTF="${ANNOTATION_DIR}/PA14_ZC01.gtf"

# Input BAM files
BAM_DIR="${PROJECT_DIR}/miseq_data/3.mapping/ZC01/BAM"

# Output directories
COUNT_DIR="${PROJECT_DIR}/miseq_data/4.gene_counting/ZC01"
RAW_COUNT_DIR="${COUNT_DIR}/raw_counts"
SUMMARY_DIR="${COUNT_DIR}/summaries"
LOG_DIR="${COUNT_DIR}/logs"

mkdir -p \
    "$ANNOTATION_DIR" \
    "$RAW_COUNT_DIR" \
    "$SUMMARY_DIR" \
    "$LOG_DIR"

log() {
    echo "[$(date '+%F %T')] $*"
}

# ------------------------
# Step 1 - Combine GFF annotations
# ------------------------
log "Combining host and phage annotations..."

cat "$HOST_GFF" "$PHAGE_GFF" \
    > "$COMBINED_GFF"

# ------------------------
# Step 2 - Convert GFF to GTF
# ------------------------
log "Converting GFF to GTF..."

gffread \
    -E \
    "$COMBINED_GFF" \
    -T \
    -o "$COMBINED_GTF" \
    > "${LOG_DIR}/gffread.log" 2>&1

# ------------------------
# Step 3 - Run featureCounts
# ------------------------
log "Starting gene counting..."

for BAM in "${BAM_DIR}"/*.bam; do

    sample=$(basename "$BAM" .bam)

    log "Counting reads for sample: ${sample}"

    featureCounts \
        -p \
        --countReadPairs \
        -T "$THREADS" \
        -s "$STRANDNESS" \
        -a "$COMBINED_GTF" \
        -t "$FEATURE_TYPE" \
        -g "$ATTRIBUTE_TYPE" \
        -o "${RAW_COUNT_DIR}/${sample}.featureCounts.tsv" \
        "$BAM" \
        > "${LOG_DIR}/${sample}.featureCounts.log" 2>&1

done

# ------------------------
# Step 4 - Organize summary files
# ------------------------
log "Collecting featureCounts summary files..."

mv "${RAW_COUNT_DIR}"/*.summary \
   "$SUMMARY_DIR"/ 2>/dev/null || true

log "Gene counting completed successfully."
