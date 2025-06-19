#!/bin/bash
#
# SLURM batch script for RNA-Seq quantification using Kallisto.
# This script performs pseudo-quantification of RNA-Seq samples using Kallisto.
# Author: Vasiliki Mantzana Oikonomaki
# License: GNU
# Last Updated: 2025-04-04

# -------------------- CONFIGURABLE VARIABLES --------------------

# Set directories (can be set as environment variables)
WORKDIR="/scratch-emmy/projects/nib00033/vicentei/quantification"
TRANSCRIPTOME="/scratch-emmy/projects/nib00033/vicentei/assemblies/EvidentialGene/output/transcript_collection.okay.mrna"
READSDIR="/scratch-emmy/projects/nib00033/vicentei/trimmed"

# -------------------- CREATE NECESSARY DIRECTORIES --------------------

mkdir -p "$WORKDIR/kallisto_quants"
cd "$WORKDIR/kallisto_quants"

# -------------------- BUILD Kallisto INDEX --------------------

# Build the Kallisto index for the transcriptome
kallisto index -i "$WORKDIR/vicentei.idx" "$TRANSCRIPTOME"

# -------------------- LIST SAMPLES --------------------

# List all sample names by extracting the base names from the fastq files
samples=$(ls $READSDIR/*_R1_001.trim.fastq.gz | sed "s/_R1_001.trim.fastq.gz//g" | uniq)

# -------------------- CREATE SAMPLE DIRECTORIES --------------------

# Make directories for each sample to store quantification results
for sample in $samples; do
    mkdir -p "$WORKDIR/kallisto_quants/$(basename $sample)"
done

# -------------------- PSEUDO-QUANTIFICATION WITH Kallisto --------------------

# Run Kallisto quantification for all samples in parallel
cd "$WORKDIR"
mkdir -p tmp

parallel --tmpdir "$WORKDIR/tmp" -j 24 kallisto quant -t 24 -i "$WORKDIR/vicentei.idx" \
    -o kallisto_quants/$(basename {}) --plaintext -b 100 "$READSDIR/{}_R1_001.trim.fastq.gz" "$READSDIR/{}_R2_001.trim.fastq.gz" ::: $samples

# -------------------- RENAME FOLDERS TO REMOVE SAMPLE NUMBERS --------------------

# List of problematic sample names to rename
declare -A rename_map=(
    ["213-CEI1L"]="CEI1L"
    ["221-CEI5L"]="CEI5L"
    ["222-CEI5S"]="CEI5S"
    ["224-CEI2L"]="CEI2L"
    ["231-CEI3L"]="CEI3L"
    ["233-CEI7S"]="CEI7S"
    ["234-CEI7L"]="CEI7L"
    ["237-CEI4L"]="CEI4L"
    ["239-CEI8L"]="CEI8L"
    ["240-CEI8S"]="CEI8S"
    ["243-EMP1L"]="EMP1L"
    ["246-EMP2L"]="EMP2L"
    ["249-EMP3L"]="EMP3L"
    ["252-EMP4L"]="EMP4L"
    ["254-EMP5S"]="EMP5S"
    ["255-EMP5L"]="EMP5L"
    ["256-LOM1S"]="LOM1S"
    ["257-LOM1L"]="LOM1L"
    ["258-LOM1E"]="LOM1E"
    ["259-LOM2S"]="LOM2S"
    ["260-LOM2L"]="LOM2L"
    ["261-LOM2E"]="LOM2E"
    ["264-LOM3E"]="LOM3E"
    ["266-CAL1L"]="CAL1L"
    ["267-CAL1E"]="CAL1E"
    ["268-CAL2S"]="CAL2S"
    ["269-CAL2L"]="CAL2L"
    ["270-CAL2E"]="CAL2E"
    ["272-CAL3L"]="CAL3L"
    ["273-CAL3E"]="CAL3E"
    ["274-CAL4S"]="CAL4S"
    ["275-CAL4L"]="CAL4L"
    ["276-CAL4E"]="CAL4E"
    ["277-LOM4S"]="LOM4S"
    ["278-LOM4L"]="LOM4L"
    ["279-LOM4E"]="LOM4E"
    ["ARO02_212-CEI1S"]="CEI1S"
    ["ARO02_225-CEI2S"]="CEI2S"
    ["ARO02_230-CEI3S"]="CEI3S"
    ["ARO02_236-CEI4S"]="CEI4S"
    ["ARO02_242-EMP1S"]="EMP1S"
    ["ARO02_245-EMP2S"]="EMP2S"
    ["ARO02_248-EMP3S"]="EMP3S"
    ["ARO02_251-EMP4S"]="EMP4S"
)

# Rename the folders based on the rename_map
for old_name in "${!rename_map[@]}"; do
    mv "$WORKDIR/kallisto_quants/$old_name" "$WORKDIR/kallisto_quants/${rename_map[$old_name]}"
done

# -------------------- CLEAN UP --------------------

# Remove temporary files
rm -rf "$WORKDIR/tmp"

# -------------------- NOTIFY COMPLETION --------------------

echo "Kallisto quantification completed successfully!" >> "$WORKDIR/kallisto_quantification_summary.txt"
