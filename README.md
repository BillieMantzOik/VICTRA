# RNA-Seq Transcript Processing Pipeline
**Overview:**
This repository contains a series of bash scripts designed to process RNA-Seq transcript data for *Oophaga vicentei*. The scripts perform transcript assembly, quantification, and post-processing using tools like SOAPdenovo, OASES, SPADES, Trinity, Idba, Kallisto, and EvidentialGene.

## Introduction
The goal of this pipeline is to take raw RNA-Seq data (FASTQ files), assemble transcripts from multiple methods (SOAPdenovo, OASES, SPADES, Trinity, Idba), and run downstream analyses such as transcript quantification (using Kallisto) and annotation using EvidentialGene.

**Key Features:**
Pre-processing of raw RNA-reads using FASTP and rcorrector.
Kraken2 classification of RNA-reads.
Multi-assembler processing for *de-novo* transcriptome assembly using SOAPdenovo, OASES, SPADES, Trinity, and Idba.
Building consensus transcsriptome collection using EvidentialGene.
Transcript quantification with Kallisto.

**System Requirements:**
This pipeline was developed to run on high-performance computing (HPC) environments using SLURM job scheduler. The following are the primary software dependencies:
SLURM: For job scheduling.
Bash: For scripting.

## **Steps and softwares used**
## 1. Preprocessing
**1. preprocess.sh**
Software Requirements:
fastp: A fast all-in-one preprocessing tool for FASTQ files. Used for adapter trimming and quality filtering.
Rcorrector: A k-mer-based error correction tool for Illumina RNA-seq reads.
perl: Required to run the run_rcorrector.pl script.
Sample Sheet (samplesheet.txt): A tab-delimited file with the following columns:
Sample name  Path to Read 1 (R1)  Path to Read 2 (R2)

Example samplesheet.txt:
`sample1	/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz sample2	/path/to/sample2_R1.fastq.gz	/path/to/sample2_R2.fastq.gz`
**Description:**
This script performs the preprocessing of RNA-seq data. It trims adapters, removes low-quality reads, and corrects sequencing errors using fastp and Rcorrector.

## 2. Kraken2 classification
**2_kraken2_classification.sh**
**Description:**
This script performs taxonomic classification of RNA-seq reads using Kraken2. It classifies the sequence data into taxonomic categories based on Kraken2's reference databases.
Software Requirements:
Kraken2: A taxonomic classification tool for sequence reads.
Krakendb: A Kraken2 database containing the required reference data for taxonomic classification.
Sample Sheet (samplesheet.txt): The sample sheet file is in the same format as the previous script.


## 3. Concatennation and normalization
**3_concatennation_and_normalization.sh**
**Description:**
This script concatenates RNA-seq paired-end files and normalizes the reads using BBNorm. It is crucial for ensuring that read counts across samples are comparable.
BBNorm: Used for normalization of reads. bbmap is used in the script.
java: bbnorm requires a Java runtime.
samplesheet:
Example
`CAL	CAL_R1.unclassified.R_1.trim.cor.fq.gz	CAL_R2.unclassified.R_2.trim.cor.fq.gz
EMP	EMP_R1.unclassified.R_1.trim.cor.fq.gz	EMP_R2.unclassified.R_2.trim.cor.fq.gz
CEI	CEI_R1.unclassified.R_1.trim.cor.fq.gz	CEI_R2.unclassified.R_2.trim.cor.fq.gz
LOM	LOM_R1.unclassified.R_1.trim.cor.fq.gz	LOM_R2.unclassified.R_2.trim.cor.fq.gz
`

## 4. rRNA removal
**4_rrna_removal.sh**
**Description:**
This script removes ribosomal RNA (rRNA) contamination from RNA-seq data. SortMeRNA is used to filter out rRNA sequences, ensuring that only the non-rRNA portions of the transcriptome are analyzed.
SortMeRNA: Tool for removing rRNA contamination from RNA-seq data.

## 5-9. Transcriptome assembly via 5 assemblers

**5_trassembly_trinity.sh**
**Description:**
This script assembles RNA-seq data into transcripts using Trinity. It is used for assembling transcriptomes when a reference genome is not available.

Trinity:A de novo transcriptome assembly tool.

**6_trassembly_rnaspades.sh**
**Description:**
This script assembles RNA-seq data into transcripts using RNAspades.

python: required to run rnaspades.py
rnaspades: Trnascriptome assembly

**7_trassembly_soap.sh**
**Description:**
This script assembles RNA-seq data into transcripts using SOAP-denovo.
SOAPdenovo: Transcriptome assembly.

**8_trassembly_velvethoases.sh**
**Description:**
This script assembles RNA-seq data into transcripts using velveth-oases.
velveth-oases: transcriptome assembly.

**9_trassembly_idba.sh**
**Description:**
This script assembles RNA-seq data into transcripts using IDBA.
idba: transcriptome assembly.

## 6. Consensus assembly using evigene
**10_consensusaseembly_EG.sh**
**Description:**
This script performs consensus assembly and gene assembly using EvidentialGene. It generates a more refined and accurate set of gene models from RNA-seq data.
SeqKit: toolkit for FASTA/Q file manipulation.
EvidentialGene: evigene software for gene assembly.

## 7. Quantification
**11_quantification.sh**
**Description:**
This script quantifies gene expression using Kallisto. It processes RNA-seq data and outputs gene expression levels for downstream analysis.
Kallisto: For transcript quantification.

## 8. Limma analysis for DGE
**12_limma_analysis.R**
**Description:**
This R script performs differential gene expression (DGE) analysis using the limma package. It processes transcript quantification results from Kallisto, annotated transcripts fro entap, and identifies differentially expressed genes.

R: Required to run the limma analysis script.
tidyverse, tximport, readr, tximportData, dplyr, edgeR, limma, DESeq2: Required R libraries for differential gene expression analysis.
transcript_to_gene_map.csv: A file containing Entap annotation results mapping transcripts to genes.
abundance.tsv: A file containing Kallisto quantification results per sample.

## 9. WGCNA Analysis for Skin and Liver Samples of *O. vicentei*
**13_WGCNA_skin.R and 14_WGCN_liver.R**
These scripts perform a Weighted Gene Co-expression Network Analysis (WGCNA) on skin (13_WGCNA_skin.R) and liver (14_WGCNA_liver.R) transcriptomic data for *O. vicentei*. The analysis identifies co-expression modules, correlates them with relevant traits, and identifies hub genes for further investigation. The results are exported for downstream analyses and visualization in tools like Cytoscape.

R: Required to run WGCNA analysis
Required R packages:
`WGCNA`, `dplyr`, `tidyverse`, `ggplot2`, `knitr`
