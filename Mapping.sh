#Alignment of processed fastq files to combined reference genome from phage and host.

#combine two reference genome (Phage-Host)
cat pa14.fasta zc03.fasta > chimeric.fasta

#Create a index
bowtie2-build /home/gwlcardoso/RNAseq_ZC03_miseq/ref_genome/chimeric/chimeric.fasta /home/gwlcardoso/RNAseq_ZC03_miseq/4.mapping/index_chimeric

#mapping with Bowtie2
for R1 in /home/gwlcardoso/RNAseq_ZC03_miseq/3.filtered_fastq/trimgalore/*_R1.fq.gz; do
R2=${R1/_R1.fq.gz/_R2.fq.gz}
sample=$(basename "$R1" _R1.fq.gz)
echo "Processing sample: $sample"
bowtie2 -x /home/gwlcardoso/RNAseq_ZC03_miseq/4.mapping/index_chimeric/index_chimeric \\
        -1 "$R1" \\
        -2 "$R2" \\
        --very-sensitive \\
        --un-conc /home/gwlcardoso/RNAseq_ZC03_miseq/4.mapping/unmapped_${sample}_%.fq \\
        -p 6 \\
        -S /home/gwlcardoso/RNAseq_ZC03_miseq/4.mapping/map_chimeric_${sample}.sam
done
