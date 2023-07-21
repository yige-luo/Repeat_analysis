if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
}
library(tidyverse)

# update the repeat class (general than repeat family)

updateBED <- function(bedFile, tsvFile){
  # Combine BED files
  bed <- read.table(bedFile, header = FALSE, col.names = c("seq_name", "start", "end", "repeat_family", "score", "strand"))
  tsv <- read_tsv(tsvFile, col_names = c("repeat_family", "repeat_class", "repeat_subclass"))
  
  bed <- 
    full_join(bed, tsv, by = "repeat_family") %>% 
    select(c("seq_name", "start", "end", "repeat_class", "score", "strand"))
  
  return(bed)
}