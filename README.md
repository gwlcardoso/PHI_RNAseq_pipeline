# PHI_RNAseq_scripts

## Overview

This repository contains a curated collection of scripts developed and used in a study of time-resolved dual RNA-seq analysis of bacteriophage–host interactions. Rather than a fully automated pipeline, this repository provides modular scripts that support each stage of transcriptomic data processing and analysis.
These scripts were designed to investigate gene expression dynamics during phage infection, enabling the characterization of both host responses and phage-driven transcriptional reprogramming.

## Scope and Philosophy

This repository is intended to:
- Provide transparent access to the analytical procedures used in the study
- Enable reproducibility of key results
- Allow flexible adaptation to different experimental designs

The workflow is not fully automated, and some steps require manual execution and parameter tuning depending on the dataset.

## Key features

- Dual RNA-seq processing for phage–host systems
- Time-series analysis of infection dynamics
- Supported for paired-end Illumina data
- Integration of quality control, alignment, and statistical analysis
- Reproducible

## Analysis Structure

1. Quality Control: Raw read assessment (FastQC and Summary reports with MultiQC)
2. Preprocessing: Adapter and low-quality trimming (Cutadapt and Fastp)
3. Alignment: Mapping reads to concatenated host and phage reference genomes (Bowtie2)
4. Post-processing: SAM/BAM manipulation (Samtools)
5. Quantification: Gene-level read counting (featureCounts)
6. Temporal Expression Analysis of Phage Genes: Identification and characterization of phage gene expression dynamics across infection time points (custom Python scripts)
7. Differential Expression Analysis: Statistical analysis (DESeq2)
8. Functional Analysis: Pathway and enrichment analysis ( GO, KEGG - ClusterProfiler)

## Usage

There is no single command to run the full analysis. Instead, scripts should be executed step-by-step following the logical order of the workflow.

## Author

Guilherme Wenceslau

MSc candidate in Bioinformatics – University of São Paulo (USP)
