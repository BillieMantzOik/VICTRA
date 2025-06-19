#!/bin/bash
#
# SLURM batch script for transcriptome assembly using SOAPdenovo-Trans for different localities.
# Designed for HPC environments with SLURM workload manager.
# Author: Vasiliki Mantzana Oikonomaki
# License: GNU
# Last Updated: 2025-04-04

# -------------------- ENVIRONMENT SETUP --------------------
module load impi/2019.5
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
source ~/.bashrc

# -------------------- CONFIGURABLE VARIABLES --------------------

# Define directories (can be set as environment variables)
READS_DIR="${READS_DIR:-/path/to/preprocessed/reads}"   # Adjust this path as needed
SAMPLESHEET="${SAMPLESHEET:-/path/to/sample_sheet/vicentei_assembly.txt}"
WORKDIR="${WORKDIR:-/path/to/assembly/output}"
TEMPDIR="$WORKDIR/soap_temp"

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
        SOAPDIR="$TEMPDIR/k$kmer"
        mkdir -p "$SOAPDIR"
        cd "$SOAPDIR"

        # Create SOAP configuration file for each k-mer size
        echo -e "#maximal read length\nmax_rd_len=150\n[LIB]\n#maximal \
read length in this lib\nrd_len_cutoff=150\n#average insert size\navg_ins=230\n#\
if sequence needs to be reversed\nreverse_seq=0\n#in which part(s) the reads \
are used\nasm_flags=3\n#minimum aligned length to contigs for a reliable read \
location (at least 32 for short insert size)\nmap_len=32\n#fastq file for \
read 1\nq1=$READ1\n#fastq file for read 2\nq2=$READ2\n" > SOAP_Conf.txt

        # Run SOAPdenovo-Trans assembly with the specified k-mer value
        SOAPdenovo-Trans-127mer all -s SOAP_Conf.txt -K $kmer -o "$SOAPDIR/SOAP_k$kmer.$loc" -F -p 30

        # Copy the resulting transcripts to the main directory
        cp "$SOAPDIR/*.scafSeq" "$WORKDIR/SOAP_k$kmer.$loc.transcripts.fa"
    done
done

# Notify completion
echo "SOAP assembly complete for all localities and k-mer sizes." >> "$WORKDIR/soap_assembly_summary.txt"
