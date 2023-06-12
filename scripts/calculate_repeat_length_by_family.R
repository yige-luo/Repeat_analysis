#1. Install and load the required packages:
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  }
if (!requireNamespace("GenomicRanges", quietly = TRUE)){
  BiocManager::install("GenomicRanges")
  }
if (!requireNamespace("readr", quietly = TRUE)){
  install.packages("readr")
  }
if (!requireNamespace("dplyr", quietly = TRUE)){
  install.packages("dplyr")
  }

library(GenomicRanges)
library(readr)
library(dplyr)

# Read BED file
args <- commandArgs(trailingOnly = TRUE)

input_bed <- args[1] # e.g. "RM_49245.SunApr91814552023/families-classified.bed"
bed_data <- read.table(input_bed, header = FALSE, col.names = c("seq_name", "start", "end", "repeat_family", "score", "strand"))

# Read classified repeat family file
input_tsv <- args[2] # e.g. "RM_49245.SunApr91814552023/families-classified.tsv"
tsv_data <- read_tsv(input_tsv, col_names = c("repeat_family", "repeat_class", "repeat_subclass"))

# update the bed_data with repeat class
bed_data <- full_join(bed_data, tsv_data, by = "repeat_family")

# Create a GRanges object from the BED file data
gr <- GRanges(
  seqnames = bed_data$seq_name,
  ranges = IRanges(start = bed_data$start + 1, end = bed_data$end),
  repeat_class = bed_data$repeat_class,
  strand = bed_data$strand
)


# Get the list of unique repeat classes
repeat_classes <- unique(tsv_data$repeat_class)

# Loop through each repeat class and calculate the total non-overlapped repeat length
for (repeat_class in repeat_classes) {
  gr_subset <- gr[gr$repeat_class == repeat_class]
  merged_gr_subset <- reduce(gr_subset)
  total_length <- sum(width(merged_gr_subset))
  cat("Total non-overlapped repeat length for", repeat_class, ":", total_length, "\n")
}
