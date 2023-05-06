#!/bin/bash

# Check if the correct number of arguments is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file="$1"
output_file="$2"

awk 'BEGIN {OFS="\t"; print "##gff-version 3"} \
/^Sequence:/ {seqid = $2} \
/^[0-9]/ {attributes = "period_size="$3";copy_number="$4";consensus_size="$5";percent_matches="$6";percent_indels="$7";score="$8; print seqid, "TRF", "tandem_repeat", $1, $2, $8, ".", ".", attributes}' "$input_file" > "$output_file"
