#!/bin/bash

#SBATCH --job-name=RepeatModeler
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-user=ygeluo@ucdavis.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load any required modules or environment settings here, if necessary
module load repeatmodeler/2.0

# Check the MAFFT_BINARIES environmental variable by
$ echo $MAFFT_BINARIES

# Unset MAFFT_BINARIES
unset MAFFT_BINARIES

# name of your database
DB_NAME=prolongata_renamed_genome_db

# Run RepeatModeler
RepeatModeler -database ${DB_NAME} -engine ncbi -pa 16 -LTRStruct

