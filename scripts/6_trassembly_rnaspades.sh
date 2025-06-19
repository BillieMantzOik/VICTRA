#!/bin/bash
#
# SLURM batch script for transcriptome assembly using SPAdes for different localities.
# Designed for HPC environments with SLURM workload manager.
# Author: Vasiliki Mantzana Oikonomaki
# License: GNU
# Last Updated: 2024-04-04

# -------------------- ENVIRONMENT SETUP --------------------
module load python
# -------------------- CONFIGURABLE VARIABLES --------------------

# Define directories (can be set as environment variables)
READS_DIR="${READS_DIR:-/path/to/preprocessed/reads}" # Adjust this path as needed
SAMPLESHEET="${SAMPLESHEET:-/path/to/sample_sheet/vicentei_assembly.txt}"
WORKDIR="${WORKDIR:-/path/to/assembly/output}"
SPADESDIR="$WORKDIR/spades_temp"
OUTDIR="$WORKDIR/spades"

# -------------------- SETUP OUTPUT FOLDERS --------------------
mkdir -p "$OUTDIR" "$SPADESDIR"

# -------------------- ASSEMBLY PARAMETERS --------------------

# Locality names (modify if necessary)
localities=(LOM EMP CEI CAL)
loc="${localities[$SLURM_ARRAY_TASK_ID - 1]}"

# -------------------- PROCESSING EACH LOCALITY --------------------

for loc in $localities; do
    # Set the input read files for each locality
    READ1="$READS_DIR/$loc.trim.cor.norm.fq.gz.rRNAfiltered_fwd.fq"
    READ2="$READS_DIR/$loc.trim.cor.norm.fq.gz.rRNAfiltered_rev.fq"

    # Create a directory for the current locality
    mkdir -p "$OUTDIR/$loc"

    # Run SPAdes assembly for the current locality
    rnaspades.py -k auto -1 "$READ1" -2 "$READ2" -o "$OUTDIR/$loc" --ss rf -m 200 -t $SLURM_CPUS_PER_TASK \
        --disable-gzip-output --checkpoints "last"

    # Copy the output files to the main directory
    cp "$OUTDIR/$loc/transcripts.fasta" "$WORKDIR/spades_$loc.transcripts.fasta"
    cp "$OUTDIR/$loc/spades.log" "$WORKDIR/spades_$loc.log"
done

# Notify completion
echo "SPAdes assembly complete for all localities." >> "$WORKDIR/assembly_summary.txt"
