# ========================
# Phage-host RNA-seq read quality control and trimming pipeline
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
MINLEN=36
QUAL=20
TRIM_R2_FRONT=16
CLEAN_INTERMEDIATES=false   # change to true if you want to delete intermediate files at the end

PROJECT_DIR="/home/gwlcardoso/Projects/phi_rnaseq/nextseq_data"
RAW_DIR="${PROJECT_DIR}/data/ZC03"
OUT_DIR="${PROJECT_DIR}/2.filtered/ZC03"

# Input QC
RAW_QC_DIR="${OUT_DIR}/qc/raw_fastqc"
RAW_MQC_DIR="${OUT_DIR}/qc/raw_multiqc"

# Step 1: cutadapt
CUTADAPT_DIR="${OUT_DIR}/01_cutadapt"

# Step 2: fastp single-end trimming
FASTP_R1_DIR="${OUT_DIR}/02_fastp_R1"
FASTP_R2_DIR="${OUT_DIR}/03_fastp_R2"
FASTP_R2_TRIM_DIR="${OUT_DIR}/04_fastp_R2_trimmed"

# Step 3: repair paired reads
REPAIR_DIR="${OUT_DIR}/05_repaired"

# Step 4: final paired-end fastp
PE_FASTP_DIR="${OUT_DIR}/06_fastp_pe"

# Final QC
FINAL_QC_DIR="${OUT_DIR}/qc/final_fastqc"
FINAL_MQC_DIR="${OUT_DIR}/qc/final_multiqc"

# Logs
LOG_DIR="${OUT_DIR}/logs"

# Adapter / primer sequences
R1_INDEX_ADAPTER="CCCTGCNNNNNNNNNNAGATCGGAAGAGC"
FASTP_ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
FASTP_ADAPTER_R2="AGATCGGAAGAGC"

mkdir -p \
  "$RAW_QC_DIR" "$RAW_MQC_DIR" \
  "$CUTADAPT_DIR" \
  "$FASTP_R1_DIR" "$FASTP_R2_DIR" "$FASTP_R2_TRIM_DIR" \
  "$REPAIR_DIR" "$PE_FASTP_DIR" \
  "$FINAL_QC_DIR" "$FINAL_MQC_DIR" \
  "$LOG_DIR"

log() {
  echo "[$(date '+%F %T')] $*"
}

# ------------------------
# Step 0 - Quality control of raw reads
# ------------------------
log "Running FastQC on raw reads..."
fastqc "${RAW_DIR}"/*.fastq.gz -o "$RAW_QC_DIR" -t "$THREADS" > "${LOG_DIR}/fastqc_raw.log" 2>&1

log "Running MultiQC on raw FastQC reports..."
multiqc "$RAW_QC_DIR" -o "$RAW_MQC_DIR" > "${LOG_DIR}/multiqc_raw.log" 2>&1

# ------------------------
# Step 1 - Remove index/adaptor from R1 with cutadapt
# ------------------------
log "Starting cutadapt trimming of R1 reads..."
for R1 in "${RAW_DIR}"/*_R1_001.fastq.gz; do
  sample=$(basename "$R1" _R1_001.fastq.gz)

  log "Cutadapt: ${sample}"
  cutadapt \
    -a "$R1_INDEX_ADAPTER" \
    --cores="$THREADS" \
    -q "$QUAL" \
    -m "$MINLEN" \
    -o "${CUTADAPT_DIR}/${sample}_R1_cutadapt.fastq.gz" \
    "$R1" \
    > "${LOG_DIR}/${sample}.cutadapt.log" 2>&1
done

# ------------------------
# Step 2 - Trim R1 with fastp
# ------------------------
log "Starting fastp trimming of R1 reads..."
for R1 in "${CUTADAPT_DIR}"/*_R1_cutadapt.fastq.gz; do
  sample=$(basename "$R1" _R1_cutadapt.fastq.gz)

  log "fastp R1: ${sample}"
  fastp \
    -i "$R1" \
    -o "${FASTP_R1_DIR}/${sample}_R1_fastp.fastq.gz" \
    -q "$QUAL" \
    -l "$MINLEN" \
    --adapter_sequence="$FASTP_ADAPTER_R1" \
    --trim_poly_g \
    --thread "$THREADS" \
    --html "${FASTP_R1_DIR}/${sample}_R1_fastp.html" \
    --json "${FASTP_R1_DIR}/${sample}_R1_fastp.json" \
    > "${LOG_DIR}/${sample}.fastp_R1.log" 2>&1
done

# ------------------------
# Step 3 - Trim R2 with fastp
# ------------------------
log "Starting fastp trimming of R2 reads..."
for R2 in "${RAW_DIR}"/*_R2_001.fastq.gz; do
  sample=$(basename "$R2" _R2_001.fastq.gz)

  log "fastp R2: ${sample}"
  fastp \
    -i "$R2" \
    -o "${FASTP_R2_DIR}/${sample}_R2_fastp.fastq.gz" \
    -q "$QUAL" \
    -l "$MINLEN" \
    --adapter_sequence="$FASTP_ADAPTER_R2" \
    --trim_poly_g \
    --thread "$THREADS" \
    --html "${FASTP_R2_DIR}/${sample}_R2_fastp.html" \
    --json "${FASTP_R2_DIR}/${sample}_R2_fastp.json" \
    > "${LOG_DIR}/${sample}.fastp_R2.log" 2>&1
done

# ------------------------
# Step 4 - Trim first 16 nt from R2
# ------------------------
log "Trimming first ${TRIM_R2_FRONT} nt from R2..."
for R2 in "${FASTP_R2_DIR}"/*_R2_fastp.fastq.gz; do
  sample=$(basename "$R2" _R2_fastp.fastq.gz)

  log "R2 front trim: ${sample}"
  fastp \
    -i "$R2" \
    -o "${FASTP_R2_TRIM_DIR}/${sample}_R2_fastp2.fastq.gz" \
    -q "$QUAL" \
    -l "$MINLEN" \
    --trim_front1 "$TRIM_R2_FRONT" \
    --trim_poly_g \
    --thread "$THREADS" \
    --html "${FASTP_R2_TRIM_DIR}/${sample}_R2_fastp2.html" \
    --json "${FASTP_R2_TRIM_DIR}/${sample}_R2_fastp2.json" \
    > "${LOG_DIR}/${sample}.fastp_R2_fronttrim.log" 2>&1
done

# ------------------------
# Step 5 - Repair paired-end reads
# ------------------------
log "Repairing paired-end reads with BBMap repair.sh..."
for R1 in "${FASTP_R1_DIR}"/*_R1_fastp.fastq.gz; do
  sample=$(basename "$R1" _R1_fastp.fastq.gz)
  R2="${FASTP_R2_TRIM_DIR}/${sample}_R2_fastp2.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "[WARN] Missing R2 for sample ${sample}: ${R2}" >&2
    continue
  fi

  log "repair.sh: ${sample}"
  repair.sh \
    in1="$R1" \
    in2="$R2" \
    out1="${REPAIR_DIR}/${sample}_R1.fixed.fastq.gz" \
    out2="${REPAIR_DIR}/${sample}_R2.fixed.fastq.gz" \
    outsingle="${REPAIR_DIR}/${sample}_singletons.fastq.gz" \
    tossbrokenreads=t \
    overwrite=t \
    > "${LOG_DIR}/${sample}.repair.log" 2>&1
done

# ------------------------
# Step 6 - Final paired-end filtering with fastp
# ------------------------
log "Running final paired-end fastp..."
for R1 in "${REPAIR_DIR}"/*_R1.fixed.fastq.gz; do
  sample=$(basename "$R1" _R1.fixed.fastq.gz)
  R2="${REPAIR_DIR}/${sample}_R2.fixed.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "[WARN] Missing repaired R2 for sample ${sample}: ${R2}" >&2
    continue
  fi

  log "fastp PE: ${sample}"
  fastp \
    -i "$R1" \
    -I "$R2" \
    -o "${PE_FASTP_DIR}/${sample}_R1_fastp.fq.gz" \
    -O "${PE_FASTP_DIR}/${sample}_R2_fastp.fq.gz" \
    --unpaired1 "${PE_FASTP_DIR}/${sample}_R1.unpaired.fq.gz" \
    --unpaired2 "${PE_FASTP_DIR}/${sample}_R2.unpaired.fq.gz" \
    -x \
    -q "$QUAL" \
    -l "$MINLEN" \
    --thread "$THREADS" \
    --html "${PE_FASTP_DIR}/${sample}_report.html" \
    --json "${PE_FASTP_DIR}/${sample}_report.json" \
    > "${LOG_DIR}/${sample}.fastp_PE.log" 2>&1
done

# ------------------------
# Step 7 - Final QC
# ------------------------
log "Running FastQC on final trimmed reads..."
fastqc "${PE_FASTP_DIR}"/*_R1_fastp.fq.gz "${PE_FASTP_DIR}"/*_R2_fastp.fq.gz \
  -o "$FINAL_QC_DIR" -t "$THREADS" > "${LOG_DIR}/fastqc_final.log" 2>&1

log "Running MultiQC on final QC reports..."
multiqc "$FINAL_QC_DIR" -o "$FINAL_MQC_DIR" > "${LOG_DIR}/multiqc_final.log" 2>&1

# ------------------------
# Optional cleanup
# ------------------------
if [[ "$CLEAN_INTERMEDIATES" == "true" ]]; then
  log "Removing intermediate files..."
  rm -f "${CUTADAPT_DIR}"/*_R1_cutadapt.fastq.gz
  rm -f "${FASTP_R1_DIR}"/*_R1_fastp.fastq.gz
  rm -f "${FASTP_R2_DIR}"/*_R2_fastp.fastq.gz
  rm -f "${FASTP_R2_TRIM_DIR}"/*_R2_fastp2.fastq.gz
fi

log "Pipeline finished successfully."
