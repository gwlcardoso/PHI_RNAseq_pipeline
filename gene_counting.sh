#combine phage-host annotations files into one
cat host.gff phage.gff > ph_combined.gff

#convert gff into gtf file
gffread -E file.gff -T -o file_name.gtf

#remove host rRNA and tRNA from .gtf file
grep -wv rRNA file.gtf > cleanrRNA_ph.gtf
grep -wv tRNA file.gtf > cleantRNA_ph.gtf

#combine to gtf cleaned files into one
cat cleanrRNA_ph.gtf  cleantRNA_ph.gtf > ph_cleaned.gtf

#quantify gene by featureCounts
featureCounts -T 6 --countReadPairs -t CDS -g transcript_id -a /path to gtf file/.gtf -o counts_ph.tsv -p -s 0 /path to bam files/*.bam
