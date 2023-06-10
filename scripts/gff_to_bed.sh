#!/bin/bash

# Check for correct number of arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_gff_file> <output_bed_file>"
    exit 1
fi

# Define input and output files
input_file="$1"
output_file="$2"

# Use awk to process the file
awk 'BEGIN {OFS="\t"} \
!/^#/ {print $1, $4, $5+1, $3, $6, $7}' "$input_file" > "$output_file"
