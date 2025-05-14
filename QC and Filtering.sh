#Commands and tools applied in step 1 - Quality Control and Filtering

#Raw data FASTQC
fastqc *.fastq.gz -o /path_to_output/ -t 6

#MultiQC
multiqc .

#Trimmomatic
for R1 in *_R1_001.fastq.gz; do
    R2="${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    sample=$(basename "$R1" _R1_001.fastq.gz)

    echo "Processing sample: $sample"

    trimmomatic PE -threads 6 \
        "$R1" "$R2" \
        "${sample}_R1_trimmed_paired.fastq.gz" "${sample}_R1_trimmed_unpaired.fastq.gz" \
        "${sample}_R2_trimmed_paired.fastq.gz" "${sample}_R2_trimmed_unpaired.fastq.gz" \
        ILLUMINACLIP:/home/gwlcardoso/RNAseq_ZC03_miseq/zc03_adapter_stub.fa:2:30:10 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50 AVGQUAL:20
        #LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 AVGQUAL:20
done

#Trimmed FASTQC
fastqc *.fastq.gz -o /path_to_output/ -t 6

#MultiQC
multiqc .

#TrimGalore
for R1 in /home/gwlcardoso/RNAseq_ZC03_miseq/1.raw_data/*_R1_001.fastq.gz; do
    # Derive the matching R2 file by replacing part of the filename
    R2=${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}
    
    # Get the sample name
    sample=$(basename "$R1" _R1_001.fastq.gz)
    
    echo "Processing sample: $sample"
    echo "R1: $R1"
    echo "R2: $R2"
    
    # Run Trim Galore
    trim_galore --paired \
        --illumina \
        --quality 20 \
        --length 50 \
        --fastqc \
        --cores 8 \
        --gzip \
        "$R1" "$R2"
done

#MultiQC
multiqc .

#Fastp
for R1 in /home/gwlcardoso/RNAseq_ZC03_miseq/1.raw_data/*_R1_001.fastq.gz; do
	R2=${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}
	sample=$(basename $R1 _R1_001.fastq.gz)
	
	fastp \
	-i "$R1" \
	-I "$R2" \
	-o "${sample}_R1_trimmed.fastq.gz" \
	-O "${sample}_R2_trimmed.fastq.gz" \
	-q 20 \
	-M 20 \
	-l 50 \
	-W 4 \
	-x \
	-5 \
	-3 \
	--adapter_fasta /home/gwlcardoso/RNAseq_ZC03_miseq/TruSeq3-PE-2.fa \
	--trim_poly_x \
	--thread 4 \
	--html "${sample}_report.html"
done

#Trimmed FASTQC
fastqc *.fastq.gz -o /path_to_output/ -t 6

#MultiQC
multiqc .
