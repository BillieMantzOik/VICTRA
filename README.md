# RNA-Seq Transcript Processing Pipeline
**Overview:**
This repository contains the code and data for the VICTRA project, focused on transcriptomic analysis of *O. vicentei* species. It includes scripts for preprocessing, assembly, classification, differential gene expression analysis and gene-co expression network analysis as well as data files (e.g., metadata and counts data) required for these analyses.

Contents:
Scripts: Preprocessing, kraken2 classification, transcript assembly, quantification, DGE analysis and WGCNA analysis

Data: Input files and results from analysis (e.g., metadata).

## Introduction
The goal of this pipeline is to take raw RNA-Seq data (FASTQ files), assemble transcripts from multiple methods (SOAPdenovo, OASES, SPADES, Trinity, Idba), and run downstream analyses such as transcript quantification (using Kallisto) and DGE analysis using limma.

**Key Features:**
Pre-processing of raw RNA-reads using FASTP and rcorrector.
Kraken2 classification of RNA-reads.
Multi-assembler processing for *de-novo* transcriptome assembly using SOAPdenovo, OASES, SPADES, Trinity, and Idba.
Building consensus transcsriptome collection using EvidentialGene.
Transcript quantification with Kallisto.
DGE analysis with limma.
Gene co-expression networks build with WGCNA.
