# VICTRA
Description
This repository contains all scripts and workflows used for:

Transcriptome assembly (de novo for O. vicentei)

Differential gene expression (DGE) analysis (limma-voom, edgeR)

Co-expression network analysis (WGCNA, Cytoscape integration)

Developed for the study Gene expression and histological changes underlying coloration differences in the Panamanian poison frog Oophaga vicentei (submitted/in press at Molecular Ecology), this pipeline ensures reproducibility of bioinformatics analyses. Raw sequencing data are available under NCBI BioProject PRJNA1234489.

Key Contents
assembly/: Scripts for quality control, assembly (e.g., Trinity/SPAdes, etc.), and assessment (BUSCO/RNAquast).

DGE/: R scripts for read quantification (Kallisto), DG statistical analysis (limma-voom), and visualization.

network_analysis/: Co-expression network construction (e.g., WGCNA) and functional enrichment (GO).

docs/: Detailed READMEs, parameter files, and example input/output.

Usage
Install dependencies (see environment.yml or requirements.txt).

Run scripts sequentially via Snakemake/Nextflow or manually (see workflow diagram).

Adjust paths/parameters in config/ files for your data.
