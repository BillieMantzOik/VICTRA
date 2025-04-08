#!/bin/bash
#
# SLURM batch script to normalize paired-end reads using BBNorm, with optional concatenation step.
# Designed for HPC environments with SLURM workload manager.
# Author: Vasiliki Mantzana Oikonomaki
# License: MIT
# Last Updated: 2025-04-04

#SBATCH --job-name=bbnorm_with_concatenation
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=10
#SBATCH --array=1-4
#SBATCH -p your_partition                        # e.g., large40
#SBATCH --output=slurm_output/slurm-norm-%x-%A_%a.out

# -------------------- ENVIRONMENT SETUP --------------------

module load gcc/9.3.0
module load openmpi/gcc.9/3.1.5
module load java/16
source ~/.bashrc

# Set OpenMP threads (adjust as needed)
export OMP_NUM_THREADS=40

# -------------------- CONFIGURABLE VARIABLES --------------------

# Define directories (can be set as environment variables)
WORKDIR="${WORKDIR:-/path/to/concatenated_reads}"
SAMPLESHEET="${SAMPLESHEET:-samplesheet_norm.txt}"

# Navigate to working directory
cd "$WORKDIR" || { echo "WORKDIR not found: $WORKDIR"; exit 1; }

# -------------------- CONCATENATION STEP --------------------

# Concatenate calovebora samples
 cat /path/to/kraken/unclassified/vicentei/*CAL*.unclassified-out.R_1.trim.cor.fq.gz \
 > /path/to/concatenated/vicentei/CAL_R1.unclassified.R_1.trim.cor.fq.gz
 cat /path/to/kraken/unclassified/vicentei/*CAL*.unclassified-out.R_2.trim.cor.fq.gz \
 > /path/to/concatenated/vicentei/CAL_R2.unclassified.R_2.trim.cor.fq.gz

# Concatenate empalizada samples
 cat /path/to/kraken/unclassified/vicentei/*EMP*.unclassified-out.R_1.trim.cor.fq.gz \
 > /path/to/concatenated/vicentei/EMP_R1.unclassified.R_1.trim.cor.fq.gz
 cat /path/to/kraken/unclassified/vicentei/*EMP*.unclassified-out.R_2.trim.cor.fq.gz \
 > /path/to/concatenated/vicentei/EMP_R2.unclassified.R_2.trim.cor.fq.gz

# Concatenate ceiba samples
 cat /path/to/kraken/unclassified/vicentei/*CEI*.unclassified-out.R_1.trim.cor.fq.gz \
 > /path/to/concatenated/vicentei/CEI_R1.unclassified.R_1.trim.cor.fq.gz
 cat /path/to/kraken/unclassified/vicentei/*CEI*.unclassified-out.R_2.trim.cor.fq.gz \
 > /path/to/concatenated/vicentei/CEI_R2.unclassified.R_2.trim.cor.fq.gz

# Concatenate loma grande samples
 cat /path/to/kraken/unclassified/vicentei/*LOM*.unclassified-out.R_1.trim.cor.fq.gz \
 > /path/to/concatenated/vicentei/LOM_R1.unclassified.R_1.trim.cor.fq.gz
 cat /path/to/kraken/unclassified/vicentei/*LOM*.unclassified-out.R_2.trim.cor.fq.gz \
 > /path/to/concatenated/vicentei/LOM_R2.unclassified.R_2.trim.cor.fq.gz

# -------------------- SAMPLE INFO EXTRACTION --------------------

# Read sample info from the sample sheet
samplename=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $1}')
f1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $2}')
f2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLESHEET" | awk '{print $3}')

echo "Running BBNorm for task ${SLURM_ARRAY_TASK_ID}: ${samplename}, files: ${f1}, ${f2}" >> normalization_array_log.txt

# -------------------- NORMALIZATION WITH BBNORM --------------------

# Find paired files automatically
f1_path="${WORKDIR}/${f1}"
f2_path="${WORKDIR}/${f2}"

bbnorm.sh -Xmx300g \
  in1="$f1_path" \
  in2="$f2_path" \
  hist="${samplename}_input.hist" \
  histout="${samplename}_output.hist" \
  out1="${samplename}_R1.trim.cor.norm.fq.gz" \
  out2="${samplename}_R2.trim.cor.norm.fq.gz" \
  target=50 \
  mindepth=2
