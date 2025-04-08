#!/bin/bash
#
# SLURM batch script to preprocess paired-end reads using fastp and Rcorrector
# Designed for HPC clusters using SLURM job scheduler.
# Designed for processing O. vicentei RNA sequence reads, provided at SRA bioproject PRJNA1234489.
# Customize paths, modules, and resource allocations as needed.
#
# Author: Vasiliki Mantzana Oikonomaki
# License: MIT
# Last Updated: 2025-04-04
# HPC cluster: NHR-NORD@Göttingen system “Emmy”

#SBATCH --job-name=preprocess_reads
#SBATCH --nodes=2
#SBATCH --ntasks=88
#SBATCH --ntasks-per-node=44
#SBATCH --cpus-per-task=4
#SBATCH --array=1-44					   	 # Job array for parallel processing of 44 samples
#SBATCH --time=12:00:00
#SBATCH -p standard                          # Replace with your HPC partition
#SBATCH --output=slurm_output/slurm-%x-%A_%a.out

# -------------------- CONFIGURABLE PARAMETERS --------------------

# Set these environment variables before running or replace with hardcoded paths
READS_DIR="${READS_DIR:-/path/to/reads}"
SAMPLESHEET="${SAMPLESHEET:-samplesheet.txt}"
OUTPUT_DIR="${OUTPUT_DIR:-./preprocessed}"
THREADS="${THREADS:-176}"

# Load required modules (adapt based on your environment)
module load openmpi/gcc.9/3.1.5
module load perl/5.22.0

# HPC environment settings
export SLURM_CPU_BIND=none
export OMP_NUM_THREADS=$THREADS

# -------------------- SCRIPT LOGIC --------------------

cd "$READS_DIR" || { echo "READS_DIR not found: $READS_DIR"; exit 1; }

# Read sample info from the samplesheet
samplename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $1}')
r1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $2}')
r2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $3}')

echo "Processing sample ${samplename}: R1=${r1}, R2=${r2}"

# Run fastp for quality filtering and adapter trimming
fastp --thread "$THREADS" --detect_adapter_for_pe \
  --in1 "$r1" \
  --in2 "$r2" \
  --out1 "${r1%%.fastq.gz}.trim.fastq.gz" \
  --out2 "${r2%%.fastq.gz}.trim.fastq.gz" \
  --html "${samplename}.fastp.html" \
  --json "${samplename}.fastp.json"

# Run Rcorrector for error correction
perl /path/to/run_rcorrector.pl \
  -t "$THREADS" \
  -1 "${r1%%.fastq.gz}.trim.fastq.gz" \
  -2 "${r2%%.fastq.gz}.trim.fastq.gz" \
  -verbose \
  -od "$OUTPUT_DIR"