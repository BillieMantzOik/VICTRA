#!/bin/bash
#
# SLURM batch script for transcriptome assembly using IDBA-UD for different localities.
# Designed for HPC environments with SLURM workload manager.
# Author: Vasiliki Mantzana Oikonomaki
# License: MIT
# Last Updated: 2025-04-04

# -------------------- CONFIGURABLE VARIABLES --------------------

# Define directories (can be set as environment variables)
READS_DIR="${READS_DIR:-/path/to/preprocessed/reads}"
SAMPLESHEET="${SAMPLESHEET:-/path/to/sample_sheet/vicentei_assembly.txt}"
WORKDIR="${WORKDIR:-/path/to/assembly/output}"
TEMPDIR="$WORKDIR/idba_temp"

# -------------------- CREATE NECESSARY DIRECTORIES --------------------

mkdir -p "$WORKDIR" "$TEMPDIR"

# -------------------- ASSEMBLY PARAMETERS --------------------

localities="LOM EMP CEI CAL"
threads=40

# -------------------- ASSEMBLY PROCESS --------------------

# Process each sample in the samplesheet using SLURM_ARRAY_TASK_ID
samplename=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$SAMPLESHEET" | awk '{print $1}')
READ1=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$SAMPLESHEET" | awk '{print $2}')
READ2=$(sed -n "$SLURM_ARRAY_TASK_ID"p "$SAMPLESHEET" | awk '{print $3}')

for loc in $localities; do
    # Create merged file for the paired-end reads
    MERGEDREADS="${READ1%%.trim.cor.norm.fq.gz.rRNAfiltered_fwd.fq}.merged.fq"
    /home/nibvasmo/idba/bin/fq2fa --merge --filter "$READ1" "$READ2" "$MERGEDREADS"

    # Run IDBA for transcript assembly
    ./idba_tran -o "$TEMPDIR" -l "$MERGEDREADS" --num_threads "$threads" --mink 19 --maxk 73 --step 18 --max_isoforms 50

    # Copy results to the main output directory
    cp "$TEMPDIR/transcript-73.fa" "$WORKDIR/idba_${loc}_transcripts.fa"
    cp "$TEMPDIR/log" "$WORKDIR/idba_${loc}_log.txt"
    
    # Clean up temporary directory for the next iteration
    rm -rf "$TEMPDIR"
done

# Notify completion
echo "IDBA assembly complete for all localities." >> "$WORKDIR/idba_assembly_summary.txt"
