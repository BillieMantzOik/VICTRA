#!/bin/bash
#
# SLURM batch script for RNA-Seq transcript processing with EvidentialGene.
# This script processes various transcript assemblies using EvidentialGene's tools.
# Author: Vasiliki Mantzana Oikonomaki
# License: MIT
# Last Updated: 2025-04-04

# -------------------- CONFIGURABLE VARIABLES --------------------

EVIGENE_DIR="/scratch/projects/nib00033/vicentei/assemblies/EvidentialGene"
ASSEMBLIES_DIR="/scratch/projects/nib00033/vicentei/assemblies"
MY_SPECIES="vicentei"

# -------------------- TRANSCRIPT PREPROCESSING --------------------

# Change to the EvidentialGene directory
cd $EVIGENE_DIR

# Prepare transcript files for EviGene
# Prepare SOAP-denovo transcripts
 /evigene/scripts/rnaseq/trformat.pl -pre $MY_SPECIES -format=soapt \
 -out $ASSEMBLIES_DIR/soap/transcripts/$MY_SPECIES"_SOAP_transcripts.fasta" -log \
 -in $ASSEMBLIES_DIR/soap/transcripts/"SOAP_k"*".transcrips.fa"

# Remove duplicates with seqkit
 seqkit rmdup -s $ASSEMBLIES_DIR/soap/transcripts/$MY_SPECIES"_SOAP_transcripts.fasta" > $EVIGENE_DIR/$MY_SPECIES"_SOAP_transcripts.fa"

# Prepare OASES transcripts
 /evigene/scripts/rnaseq/trformat.pl -pre $MY_SPECIES -format=velvet \
 -out $EVIGENE_DIR/$MY_SPECIES"_OASES_transcripts.fasta" -log \
 -in $ASSEMBLIES_DIR/oases/transcripts/"oases_k"*"transcrips.fasta"

# Remove duplicates with seqkit
 seqkit rmdup -s $EVIGENE_DIR/$MY_SPECIES"_OASES_transcripts.fasta" > $EVIGENE_DIR/$MY_SPECIES"_oases_transcripts.fa"

# Prepare SPADES transcripts
 /evigene/scripts/rnaseq/trformat.pl -pre vicenteispad \
 -out $EVIGENE_DIR/$MY_SPECIES"_SPADES_transcripts.fa" -log \
 -in $ASSEMBLIES_DIR/spades/"transcripts_"*".fasta"

# Prepare Trinity transcripts
 /evigene/scripts/rnaseq/trformat.pl -pre $MY_SPECIES -format=trinity \
 -out $EVIGENE_DIR/$MY_SPECIES"_Trinity_transcripts.fa" -log \
 -in $ASSEMBLIES_DIR/trinity/"trinity."*".Trinity.fasta"

# Prepare Idba_tran transcripts
 /evigene/scripts/rnaseq/trformat.pl -pre $MY_SPECIES -format=idba \
 -out $EVIGENE_DIR/$MY_SPECIES"_Idba_transcripts.fa" -log \
 -in $ASSEMBLIES_DIR/idba/transcripts/"idba_"*"transcripts.fa"

# -------------------- CONCATENATE ALL TRANSCRIPTS --------------------

# Concatenate all processed transcript files into one collection
 cat $EVIGENE_DIR/*.fa > $EVIGENE_DIR/$MY_SPECIES"_transcript_collection.fa"

# -------------------- RUN TR2AACDS --------------------

# Run the tr2aacds4.pl script to process the concatenated transcript collection
/evigene/scripts/prot/tr2aacds4.pl -tidy -NCPU 32 -MAXMEM 640880 -MINAA 70 -reorient -species Oophaga_vicentei -pHeterozygosity=2 -log -cdna $EVIGENE_DIR/vicentei_transcript_collection.fa

# -------------------- CLEAN UP --------------------

# Clean up any temporary files if necessary
# rm -f $EVIGENE_DIR/tmp/*

# -------------------- NOTIFY COMPLETION --------------------

echo "EvidentialGene processing completed successfully!" >> "$EVIGENE_DIR/vicentei_oa_summary.txt"
