#!/bin/bash

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <stk_file> <output_bed_file>"
  exit 1
fi

stk_file="$1"
output_bed_file="$2"

awk '
  /^#=GF ID/ { family = $3 }
  /^#=GF SQ/ { score = $3 }
  /^[^#]/ && NF > 0 && !/^\/\// {
    split($1, coords, "[:-]")
    if (coords[2] > coords[3]) {
      chrom_start = coords[3]
      chrom_end = coords[2]
      strand = "-"
    } else {
      chrom_start = coords[2]
      chrom_end = coords[3]
      strand = "+"
    }
    printf "%s\t%d\t%d\t%s\t%d\t%s\n", coords[1], chrom_start, chrom_end, family, score, strand
  }
' "$stk_file" > "$output_bed_file"
