#!/bin/bash
#JOB NAME
#SBATCH --job-name=krakenlb
#ACOUNT
#SBATCH -A nib00033
#RESOURCES ALLOCATION
#PARTITION
#SBATCH -p standard96

#LOG SCRIPT mail
#SBATCH --output=/home/nibvasmo/scripts/slurmOut/slurm-%x.oe
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=vasiliki.mantzana.oikonomaki@tiho-hannover.de

#ENVIRONMENT
#load environemt
export PREFERRED_SOFTWARE_STACK=nhr-lmod
source /sw/etc/profile/profile.sh
source ~/.bashrc
module load gcc
module load perl

KRAKENDB=/scratch-emmy/projects/nib00033/Vasiliki/microbiome
cd /scratch-emmy/projects/nib00033/Vasiliki/
#kraken library
kraken2-build --build --db $KRAKENDB