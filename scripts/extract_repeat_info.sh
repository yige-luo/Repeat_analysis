#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <classified_file> <output_tsv_file>"
  exit 1
fi

classified_file="$1"
output_tsv_file="$2"

# Extract relevant information and print as a tab-delimited table
grep "^>" "$classified_file" | 
awk -F '[#/_]' '/^>/ {
  repeat_family = $1 "_" $2;
  gsub(/^>/, "", repeat_family);
  repeat_class = $3;
  gsub(/ .*$/, "", repeat_class);
  repeat_subtype = $4;
  gsub(/ .*$/, "", repeat_subtype);
  printf("%s\t%s\t%s\n", repeat_family, repeat_class, repeat_subtype);
}' > "$output_tsv_file"