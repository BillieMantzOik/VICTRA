#!/bin/bash
#
# SLURM batch script for transcriptome assembly using Trinity for different localities.
# Designed for HPC environments with SLURM workload manager.
# Author: Vasiliki Mantzana Oikonomaki
# License: MIT
# Last Updated: 2024-04-04

# -------------------- ENVIRONMENT SETUP --------------------
#ENVIRONMENT
module load impi/2019.5
module load gcc/9.3.0
export OMP_NUM_THREADS=40
source ~/.bashrc
LD_LIBRARY_PATH=/home/nibvasmo/rcorrector/jellyfish/.libs:$LD_LIBRARY_PATH

# -------------------- CONFIGURABLE VARIABLES --------------------

READS_DIR="/scratch-emmy/projects/nib00033/preprocess/vicentei"
samplesheet="/home/nibvasmo/array_samplesheets/vicentei/vicentei_assembly.txt"
WORKDIR="/scratch/projects/nib00033/vicentei/assemblies"
TRINITYDIR="$WORKDIR/trinity"
mkdir -p $TRINITYDIR
cd $TRINITYDIR

# --- STEP 1: Sample Information ---
# Arraying: read the sample sheet and assign read names dynamically
samplename=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}')
READ1=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $2}')
READ2=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $3}')
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${samplename}, READ1 is ${READ1}, and READ2 is ${READ2}." >> output_array_trin1.txt

# --- STEP 2: Run Trinity Assembly ---
localities="LOM CEI EMP CAL"
for loc in $localities; do  
        # Run Trinity assembly for the current sample
        cd $TRINITYDIR
        Trinity --seqType fq --left $READ1 --right $READ2 --SS_lib_type RF --CPU 40 \
        --max_memory 200G --bflyCalculateCPU \
        --full_cleanup --min_kmer_cov 3 \
        --no_normalize_reads --no_version_check \
        --output "trinity_$loc"
        
        # Copy the resulting Trinity fasta file to the WORKDIR
        cp "trinity_$loc.Trinity.fasta" $WORKDIR
    done
done

echo "Trinity assembly for all localities is complete." >> $WORKDIR/trinity_assembly_status.txt