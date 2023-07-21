if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
}
library(tidyverse)
source("scripts/helper_func/update_bed_repeat_class.R")

species <- "Dpro"
# species <- "Drho"
# species <- "Dcar"
# species <- "Dmel"
analysis <- "RepeatModeler"
partition <- ""
if(species == "Dpro"){
  # partition <- "removed"
  # partition <- "deduped"
  partition <- "renamed"
}
inDate <- "2023-07-14"
# outDate <- Sys.Date()

inPath <- file.path("data/interim", species, analysis, partition, inDate)
bedFiles <- list.files(inPath, "*.bed", recursive = TRUE, full.names = TRUE)
tsvFiles <- list.files(inPath, "*.tsv", recursive = TRUE, full.names = TRUE)
outPath <- file.path("data/final", species, analysis, partition)

if(!dir.exists(outPath)){
  dir.create(outPath, showWarnings = FALSE, recursive = TRUE)
}

# check if they are in the same directory
all.equal(dirname(bedFiles), dirname(tsvFiles))

# merge repeatmodeler bed file across multiple runs
combined_bed <- 
  lapply(seq_along(bedFiles), function(i){
    updateBED(bedFiles[i], tsvFiles[i])
    }) %>% 
  do.call("rbind", .)

summary(combined_bed)

write.table(combined_bed, paste0(outPath, "/", "families-classified.bed"), 
            quote = FALSE, sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
