# For a GRange object with multiple `seqnames`, generate a figure facet wrapped by 
# `seqnames`. Inside each panel, the x-axis is the coordinate of genomic ranges, 
# the y-axis does not have actual meaning, and each feature record is spaced evenly 
# on the y-axis as a filled rectangle, color-coded by `fill_var`.


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# if (!requireNamespace("GenomicRanges", quietly = TRUE)){
#   BiocManager::install("GenomicRanges")
# }
# if (!requireNamespace("furrr", quietly = TRUE)) {
#   install.packages("furrr")
# }
# if (!requireNamespace("tidyverse", quietly = TRUE)) {
#   install.packages("tidyverse")
# }
library(GenomicRanges)
library(tidyverse)
library(furrr)

plot_GRanges <- function(gr, fill_var = "repeat_class"){
  # convert the GRanges object to a data frame
  gr_df <- as.data.frame(gr)
  
  gr_df <- gr_df %>%
    group_by(seqnames) %>%
    mutate(rank = row_number())
  
  fill_var <- rlang::sym(fill_var)
  
  # plot the 
  p <- ggplot(gr_df, aes(xmin = start, xmax = end, ymin = rank, ymax = rank+0.5, 
                      fill = !!fill_var)) +
    geom_rect() +
    facet_wrap(~ seqnames, scales = "free") +
    xlab("Coordinate") +
    ylab("") +
    labs(fill = "Repeat Class") + 
    scale_fill_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}
