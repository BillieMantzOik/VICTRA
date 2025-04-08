#vicentei RNA-Seq Transcript Processing Pipeline
Overview
This repository contains a series of bash scripts designed to process RNA-Seq transcript data for Oophaga vicentei. The scripts perform transcript assembly, quantification, and post-processing using tools like SOAPdenovo, OASES, SPADES, Trinity, Idba, Kallisto, and EvidentialGene.

Table of Contents
Introduction

System Requirements

Installation

Pipeline Workflow

Script Details

Running the Pipeline

Troubleshooting

Contact

Introduction
The goal of this pipeline is to take raw RNA-Seq data (FASTQ files), assemble transcripts from multiple methods (SOAPdenovo, OASES, SPADES, Trinity, Idba), and run downstream analyses such as transcript quantification (using Kallisto) and annotation using EvidentialGene.

Key Features:
Multi-assembler processing using SOAPdenovo, OASES, SPADES, Trinity, and Idba.

Transcript quantification with Kallisto.

Post-processing of the transcript data using EvidentialGene.

Parallelized computations to speed up the process.

System Requirements
This pipeline was developed to run on high-performance computing (HPC) environments using SLURM job scheduler. The following are the primary software dependencies:

SLURM: For job scheduling.

Bash: For scripting.

Modules: Used for loading necessary environments (e.g., impi/2019.5, gcc, hdf5-parallel).

SOAPdenovo: For transcript assembly.

OASES: For transcript assembly.

SPADES: For transcript assembly.

Trinity: For transcript assembly.

Idba: For transcript assembly.

Kallisto: For transcript quantification.

EvidentialGene: For transcript post-processing and annotation.

Installation
Ensure that the required software and modules are available on your system before running this pipeline. You may need to load specific modules in your environment using the module load command. Additionally, you may need to adjust paths based on where the software is installed on your system.

The scripts in this pipeline expect certain directory structures. Here’s an example setup:

bash
Copy
Edit
/scratch/emmy/projects/nib00033/preprocess/vicentei/       # Input reads directory
/scratch/projects/nib00033/vicentei/assemblies/              # Assemblies output directory
/home/nibvasmo/array_samplesheets/vicentei/vicentei_assembly.txt  # Samplesheet
Pipeline Workflow
The pipeline involves several key steps:

Transcript Assembly:

SOAPdenovo, OASES, SPADES, Trinity, and Idba are used for transcript assembly.

Each assembly method generates transcript files that are then processed and cleaned.

Transcript Quantification:

Kallisto is used to quantify the assembled transcripts based on a reference transcriptome.

Post-Processing with EvidentialGene:

The transcript files are processed using EvidentialGene's tr2aacds4.pl to clean, filter, and annotate them.

Script Details
Main Scripts
vicentei-soap.sh:

Runs SOAPdenovo for transcript assembly.

Parallelizes jobs based on different k-mer values.

Processes multiple locations (e.g., LOM, EMP, CEI, CAL).

Outputs assembled transcripts for each locality.

vicentei-oa.sh:

Runs the OASES assembler for transcript assembly.

Outputs assembled transcripts for each locality.

vicentei-idba.sh:

Runs the Idba-Tran assembler for transcript assembly.

Merges paired-end reads and runs the assembly.

Outputs assembled transcripts for each locality.

vicentei-kallisto.sh:

Runs Kallisto for transcript quantification.

Uses the reference transcriptome to generate pseudo-quantifications.

Outputs the results for each sample.

vicentei-evigene.sh:

Runs EvidentialGene's tr2aacds4.pl script for transcript post-processing.

Generates a final processed transcript collection.

Running the Pipeline
Step 1: Set Up Your Environment
Load the necessary modules: The pipeline relies on specific software modules such as impi, gcc, and hdf5-parallel. Ensure that they are available on your system.

Prepare input data:

Ensure your input RNA-Seq FASTQ files are placed in the correct directory (/scratch-emmy/projects/nib00033/preprocess/vicentei/).

Ensure your sample sheet is located at /home/nibvasmo/array_samplesheets/vicentei/vicentei_assembly.txt.

Directory structure:

Make sure the pipeline's output directories exist, or modify the script to create them.

Step 2: Submit Jobs
To run the pipeline, submit the SLURM jobs for each script:

For SOAPdenovo assembly:

bash
Copy
Edit
sbatch vicentei-soap.sh
For OASES assembly:

bash
Copy
Edit
sbatch vicentei-oa.sh
For Idba assembly:

bash
Copy
Edit
sbatch vicentei-idba.sh
For Kallisto quantification:

bash
Copy
Edit
sbatch vicentei-kallisto.sh
For EvidentialGene post-processing:

bash
Copy
Edit
sbatch vicentei-evigene.sh
Step 3: Monitor Progress
Monitor the progress of each SLURM job by checking the SLURM output files in the /home/nibvasmo/scripts/slurmOut/ directory. These files are generated with the job ID and script name in their filenames.

Example:

bash
Copy
Edit
/home/nibvasmo/scripts/slurmOut/slurm-vicentei-soap-12345.out
You will receive an email when the job starts and finishes if you have specified the correct email address.
