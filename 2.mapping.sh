# ========================
# RNA-seq mapping, coverage, and BigWig generation
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
REMOVE_SAM=true

PROJECT_DIR="/home/gwlcardoso/Projects/phi_rnaseq/nextseq_data"

# Reference genome
REF_FASTA="${PROJECT_DIR}/ref_gen/concat_rRNA/PA14_ZC03.fna"

# Input trimmed reads
READS_DIR="${PROJECT_DIR}/2.filtered/ZC03/repaired/fastp_pe_x"

# Output directories
MAP_DIR="${PROJECT_DIR}/3.mapping/rRNA_removed"
INDEX_DIR="${MAP_DIR}/bowtie2_index"

SAM_DIR="${MAP_DIR}/sam"
BAM_DIR="${MAP_DIR}/bam"
LOG_DIR="${MAP_DIR}/logs"
STATS_DIR="${MAP_DIR}/stats"
COVERAGE_DIR="${MAP_DIR}/coverage"
BIGWIG_DIR="${MAP_DIR}/bigwig"

mkdir -p \
  "$INDEX_DIR" \
  "$SAM_DIR" \
  "$BAM_DIR" \
  "$LOG_DIR" \
  "$STATS_DIR" \
  "$COVERAGE_DIR" \
  "$BIGWIG_DIR"

# Bowtie2 index basename
BT2_INDEX="${INDEX_DIR}/PA14_ZC03"

log() {
    echo "[$(date '+%F %T')] $*"
}

# ------------------------
# Step 1 - Build Bowtie2 index
# ------------------------
log "Building Bowtie2 index..."

bowtie2-build \
    "$REF_FASTA" \
    "$BT2_INDEX" \
    > "${LOG_DIR}/bowtie2-build.log" 2>&1

# ------------------------
# Step 2 - Map paired-end reads with Bowtie2
# ------------------------
log "Starting Bowtie2 mapping..."

for R1 in "${READS_DIR}"/*_R1_fastp.fq.gz; do

    sample=$(basename "$R1" _R1_fastp.fq.gz)
    R2="${READS_DIR}/${sample}_R2_fastp.fq.gz"

    if [[ ! -f "$R2" ]]; then
        echo "[WARN] Missing R2 file for sample ${sample}" >&2
        continue
    fi

    log "Mapping sample: ${sample}"

    bowtie2 \
        -x "$BT2_INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        --very-sensitive-local \
        --phred33 \
        --threads "$THREADS" \
        --un-conc-gz "${MAP_DIR}/unmapped_${sample}_%.fq.gz" \
        -S "${SAM_DIR}/${sample}.sam" \
        2> "${LOG_DIR}/${sample}_bowtie2.log"

done

# ------------------------
# Step 3 - Convert SAM to sorted BAM and index
# ------------------------
log "Converting SAM to sorted BAM..."

for SAM in "${SAM_DIR}"/*.sam; do

    sample=$(basename "$SAM" .sam)

    log "Processing BAM: ${sample}"

    samtools view \
        -@ "$THREADS" \
        -bS "$SAM" | \
    samtools sort \
        -@ "$THREADS" \
        -o "${BAM_DIR}/${sample}.sorted.bam"

    samtools index \
        -@ "$THREADS" \
        "${BAM_DIR}/${sample}.sorted.bam"

done

# ------------------------
# Step 4 - Generate alignment statistics
# ------------------------
log "Generating mapping statistics..."

for BAM in "${BAM_DIR}"/*.sorted.bam; do

    sample=$(basename "$BAM" .sorted.bam)

    # Flagstat
    samtools flagstat "$BAM" \
        > "${STATS_DIR}/${sample}_flagstat.txt"

    # idxstats
    samtools idxstats "$BAM" \
        > "${STATS_DIR}/${sample}_idxstats.txt"

    # stats
    samtools stats "$BAM" \
        > "${STATS_DIR}/${sample}_samtools_stats.txt"

done

# ------------------------
# Step 5 - Calculate genome coverage
# ------------------------
log "Calculating coverage..."

for BAM in "${BAM_DIR}"/*.sorted.bam; do

    sample=$(basename "$BAM" .sorted.bam)

    samtools coverage "$BAM" \
        > "${COVERAGE_DIR}/${sample}_coverage.tsv"

done

# ------------------------
# Step 6 - Generate normalized BigWig files
# ------------------------
log "Generating BigWig coverage tracks..."

for BAM in "${BAM_DIR}"/*.sorted.bam; do

    sample=$(basename "$BAM" .sorted.bam)

    # Forward strand
    bamCoverage \
        -b "$BAM" \
        -o "${BIGWIG_DIR}/${sample}_forward.bw" \
        --filterRNAstrand forward \
        --binSize 1 \
        --normalizeUsing RPKM \
        --numberOfProcessors "$THREADS"

    # Reverse strand
    bamCoverage \
        -b "$BAM" \
        -o "${BIGWIG_DIR}/${sample}_reverse.bw" \
        --filterRNAstrand reverse \
        --binSize 1 \
        --normalizeUsing RPKM \
        --numberOfProcessors "$THREADS"

done

# ------------------------
# Step 7 - Remove SAM files (optional)
# ------------------------
if [[ "$REMOVE_SAM" == "true" ]]; then

    log "Removing intermediate SAM files..."
    rm -f "${SAM_DIR}"/*.sam

fi

log "Pipeline completed successfully."
