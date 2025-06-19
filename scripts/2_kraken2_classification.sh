#!/bin/bash
#
# SLURM batch script to run Kraken2 classification on paired-end reads
# Designed for use on HPC clusters using SLURM job scheduler.
# Designed for preprocessed O. vicentei RNA sequence reads, provided at SRA bioproject PRJNA1234489.
# Author: Vasiliki Mantzana Oikonomaki
# License: GNU
# Last Updated: 2025-04-04

#SBATCH --job-name=kraken2_classification
#SBATCH --nodes=4
#SBATCH --ntasks=88
#SBATCH --ntasks-per-node=22
#SBATCH --array=1-44
#SBATCH --time=48:00:00
#SBATCH --output=slurm_output/slurm-%x-%A-%a.out


# -------------------- ENVIRONMENT SETUP --------------------

# Load required modules (adjust versions as needed)
module load gcc/9.3.0
module load openmpi/gcc.9/3.1.5

# Set OpenMP threads
export OMP_NUM_THREADS=40

# -------------------- CONFIGURABLE VARIABLES --------------------

# Input and output directories (can be passed via environment variables)
READS_DIR="${READS_DIR:-/path/to/preprocessed_reads}" # reads after preprocess
KRAKENDB="${KRAKENDB:-/path/to/kraken2_db}" # krakendb
WORKDIR="${WORKDIR:-/path/to/working_dir}"
SAMPLESHEET="${SAMPLESHEET:-samplesheet_kraken.txt}"  # txt file containing sample names and path of r1 r2 after preprocessing step.

# Move to working directory
cd "$WORKDIR" || { echo "WORKDIR not found: $WORKDIR"; exit 1; }

# -------------------- SAMPLE INFO EXTRACTION --------------------

# Read sample info from the sample sheet
samplename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $1}')
f1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $2}')
f2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $3}')

echo "Running Kraken2 for task ${SLURM_ARRAY_TASK_ID}: sample ${samplename} with reads ${f1}, ${f2}" >> kraken_array_output_log.txt

# -------------------- KRAKEN2 CLASSIFICATION --------------------

kraken2 \
  --db "$KRAKENDB" \
  --threads "$OMP_NUM_THREADS" \
  --paired \
  --classified-out "${samplename}.classified-out.R#.fq" \
  --unclassified-out "${samplename}.unclassified-out.R#.fq" \
  --gzip-compressed \
  --confidence 0.5 \
  --output "${samplename}.kraken.out" \
  --report "${samplename}.kraken.report.txt" \
  --use-names \
  "$READS_DIR/$f1" "$READS_DIR/$f2"
