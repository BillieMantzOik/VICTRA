#!/bin/bash
#
# SLURM batch script to remove rRNA contamination from paired-end reads using SortMeRNA.
# Designed for HPC environments with SLURM workload manager.
# Author: Vasiliki Mantzana Oikonomaki
# License: MIT
# Last Updated: 2024-04-04

# -------------------- ENVIRONMENT SETUP --------------------
source ~/.bashrc
# -------------------- CONFIGURABLE VARIABLES --------------------

# Define directories
READS_DIR="${READS_DIR:-/path/to/concatenated_reads}"  # Adjust this path as needed
WORKDIR="${WORKDIR:-/path/to/sortme}"
SAMPLESHEET="${SAMPLESHEET:-samplesheet_sort.txt}"  # Replace with your actual sample sheet

# Navigate to working directory
cd "$WORKDIR" || { echo "WORKDIR not found: $WORKDIR"; exit 1; }

# -------------------- SORTMERNA RRNA REMOVAL --------------------

# Iterate through the files in the READS_DIR and process each pair
for r1 in "$READS_DIR"/*_R1.trim.cor.norm.fq.gz; do
    r2="${r1%%_R1.trim.cor.norm.fq.gz}"_R2.trim.cor.norm.fq.gz

    sortmerna --ref /scratch-emmy/usr/nibvasmo/smr_v4.3_fast_db.fasta \
        --reads "$r1" \
        --reads "$r2" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --workdir "$WORKDIR" \
        --blast 1 \
        -v \
        --paired_in \
        --fastx \
        --out2 \
        --aligned "$r1.rRNAcontaminated" \
        --other "$r1.rRNAfiltered" \
        --task 2 \
        --dbg-level 1

    # Clean up temporary files
    rm -rf "$WORKDIR/kvdb"
done
