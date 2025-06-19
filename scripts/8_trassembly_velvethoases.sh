#!/bin/bash
#
# SLURM batch script for transcriptome assembly using Oases for different localities.
# Designed for HPC environments with SLURM workload manager.
# Author: Vasiliki Mantzana Oikonomaki
# License: GNU
# Last Updated: 2025-04-04

# -------------------- CONFIGURABLE VARIABLES --------------------

# Define directories (can be set as environment variables)
READS_DIR="${READS_DIR:-/path/to/preprocessed/reads}"   # Adjust this path as needed
SAMPLESHEET="${SAMPLESHEET:-/path/to/sample_sheet/vicentei_assembly.txt}"
WORKDIR="${WORKDIR:-/path/to/assembly/output}"
TEMPDIR="$WORKDIR/voases_temp"

# -------------------- SETUP OUTPUT FOLDERS --------------------
mkdir -p "$WORKDIR" "$TEMPDIR"

# -------------------- ASSEMBLY PARAMETERS --------------------

# Locality names (modify if necessary)
localities="LOM EMP CEI CAL"
kmerValues="19 37 55 73"  # List of k-mer sizes to iterate over

# -------------------- PROCESSING EACH LOCALITY --------------------

for loc in $localities; do
    # Set the input read files for each locality
    READ1="$READS_DIR/$loc.trim.cor.norm.fq.gz.rRNAfiltered_fwd.fq"
    READ2="$READS_DIR/$loc.trim.cor.norm.fq.gz.rRNAfiltered_rev.fq"

    # Iterate over k-mer values
    for kmer in $kmerValues; do
        OASESDIR="$TEMPDIR/k$kmer"
        mkdir -p "$OASESDIR"
        cd "$OASESDIR"

        # Run Velveth and Velvetg for Oases assembly
        echo "VELVETH_START"
        velveth $OASESDIR $kmer -shortPaired -fastq -separate $READ1 $READ2
        echo "VELVETG_START"
        velvetg $OASESDIR -read_trkg yes -min_contig_lgth 100 -cov_cutoff 4 -ins_length 400 -clean yes
        echo "OASES_START"
        oases $OASESDIR -cov_cutoff 4 -min_trans_lgth 200

        # Copy resulting transcripts to the main directory
        cp "$OASESDIR/transcripts.fa" "$WORKDIR/oases_k$kmer_$loc.transcripts.fa"

        # Clean up temporary directory
        rm -rf "$OASESDIR"
    done
done

# Notify completion
echo "Oases assembly complete for all localities and k-mer sizes." >> "$WORKDIR/oases_assembly_summary.txt"
