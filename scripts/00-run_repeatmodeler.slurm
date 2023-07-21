#!/bin/bash

#SBATCH --job-name=RepeatModeler
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --mail-user=ygeluo@ucdavis.edu
#SBATCH --mail-type=BEGIN,END,FAIL

if [ $# -ne 3]; then
    echo "Usage: $0 <input_FASTA_file> <repeatmodeler_database_name> <random_seed>"
    exit 1
fi
# The path of input genome assembly (FASTA_PATH)
FASTA_PATH=$1
# Define the DB_NAME
DB_NAME=$2
# Define random seed for reproducibility
SRAND=$3

# Load any required modules or environment settings here, if necessary
module load repeatmodeler/2.0

# Check the MAFFT_BINARIES environmental variable by
$ echo $MAFFT_BINARIES

# Unset MAFFT_BINARIES
unset MAFFT_BINARIES

# Check if files exist
FILES=("${DB_NAME}.nhr" "${DB_NAME}.nin" "${DB_NAME}.nnd" "${DB_NAME}.nni" "${DB_NAME}.nog" "${DB_NAME}.nsq" "${DB_NAME}.translation")
MISSING=false

for file in "${FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "$file not found!"
        MISSING=true
        break
    fi
done

# If any file is missing, run BuildDatabase
if [ "$MISSING" = true ]; then
    echo "Running BuildDatabase..."
    BuildDatabase -name ${DB_NAME} -engine ncbi ${FASTA_PATH}
else
    echo "All necessary files are present. No need to run BuildDatabase."
fi

# Run RepeatModeler
echo "Running RepeatModeler..."
RepeatModeler -database ${DB_NAME} -engine ncbi -pa 16 -srand ${SRAND} -LTRStruct


