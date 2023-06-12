# Calculate total length of tandem repeats (consider overlapping regions)

#1. Install and load the required packages:
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  }
if (!requireNamespace("rtracklayer", quietly = TRUE)){
  BiocManager::install("rtracklayer")
  }
if (!requireNamespace("GenomicRanges", quietly = TRUE)){
  BiocManager::install("GenomicRanges")
  }

library(rtracklayer)
library(GenomicRanges)

#2. Import the GFF3 file:
args <- commandArgs(trailingOnly = TRUE)

gff3_file <- args[1]
# gff3_file <- "E:/Sequencing/Novogene/C202SC19030542/2020-11-18/References/tandem_repeat_finder/Dcar/carrolli_assembly_sm_TRF.gff3"
gff3_data <- import.gff3(gff3_file)

#3. Create a GenomicRanges object from the GFF3 data:
gr <- GRanges(seqnames = seqnames(gff3_data),
              ranges = IRanges(start = start(gff3_data), end = end(gff3_data)))

#4. Reduce the GenomicRanges object to merge overlapping regions:
reduced_gr <- reduce(gr)

#5. Calculate the total length of tandem repeat regions:
total_length <- sum(width(reduced_gr))
print(sprintf("total length of tandem repeat regions: %d bp", total_length))

# 32,847,966 bp
