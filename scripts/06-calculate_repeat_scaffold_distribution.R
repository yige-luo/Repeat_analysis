# Required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!requireNamespace("Biostrings", quietly = TRUE)){
  BiocManager::install("Biostrings")
}
if (!requireNamespace("GenomicRanges", quietly = TRUE)){
  BiocManager::install("GenomicRanges")
}
if (!requireNamespace("rtracklayer", quietly = TRUE)){
  BiocManager::install("rtracklayer")
}
if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
}
if (!requireNamespace("furrr", quietly = TRUE)) {
  install.packages("furrr")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
if (!requireNamespace("progress", quietly = TRUE)) {
  install.packages("progress")
}
library(Biostrings)
library(rtracklayer)
library(GenomicRanges)
library(progress)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(furrr)
library(readr)
library(tibble)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
genomeFile <- args[1]
species <- args[2]
bedFile1 <- args[3]
bedFile2 <- args[4]

# Function to handle the requirements
processGenome <- function(genomeFile, species, bedFile1, bedFile2) {
  
  # Initialize the progress bar
  pb <- progress_bar$new(format = "[:bar] :percent eta: :eta")
  pb$tick()
  
  # Prepare for output
  analysis <- "Table_Plot_Prototype"
  partition <- ""
  if(species == "Dpro"){
    # partition <- "removed"
    # partition <- "deduped"
    partition <- "renamed"
  }
  
  msg <- sprintf("Generating %s in %s!", analysis, species)
  print(msg)
  
  outdata_dir <- file.path("data/interim", species, analysis, partition, Sys.Date())
  outresult_dir <- file.path("results", species, analysis, partition, Sys.Date())
  if(!dir.exists(outdata_dir)){
    dir.create(outdata_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if(!dir.exists(outresult_dir)){
    dir.create(outresult_dir, recursive = TRUE, showWarnings = FALSE)
  }

  print("Loading genome and genome feature files!")
  
  # Combine BED files
  bed1 <- read.table(bedFile1, 
                     col.names = c("seq_name", "start", "end", "repeat_class", "score", "strand"))
  bed2 <- read.table(bedFile2, 
                     col.names = c("seq_name", "start", "end", "repeat_class", "score", "strand"))
  combinedBed <- rbind(bed1, bed2)
  
  # Load Genome file
  genome <- readDNAStringSet(genomeFile)
  

  print("Creating GRanges objects!")
  pb$tick()
  
  # Create a GRanges object
  
  gr <- GRanges(seqnames = combinedBed$seq_name,
                ranges = IRanges(start = combinedBed$start, end = combinedBed$end - 1),
                strand = combinedBed$strand, 
                repeat_class = combinedBed$repeat_class,
                interval_width = combinedBed$end - combinedBed$start + 1)
  
  # Ignore strand information
  strand(gr) <- "*"
  
  # Separate genomic ranges by unique values of the metadata column 'repeat_class'
  gr_list <- split(gr, mcols(gr)$repeat_class)
  
  # Apply reduce to each subset of genomic ranges and assign corresponding repeat_class
  gr_reduced_list <- lapply(names(gr_list), function(repeat_class) {
    gr_reduced <- reduce(gr_list[[repeat_class]])
    mcols(gr_reduced)$repeat_class <- repeat_class
    gr_reduced
  })
  
  # Combine the reduced GRanges objects back into one GRanges object
  gr_reduced <- GRanges()
  for(g in gr_reduced_list){
    gr_reduced <- append(gr_reduced, g)
  }
  
  # Partition genomic intervals
  disjoined <- disjoin(gr_reduced)
  annotations <- findOverlaps(disjoined, gr_reduced)
  
  print("Repeat class assignment started: computational expensive!")
  print(Sys.time())
  pb$tick()
  
  # plan to use multiple cores for processing
  plan(multisession)
  
  # Assign all annotations to each partitioned interval in parallel
  disjoined_annotations <- future_map(seq_along(disjoined), function(i) {
    hits <- subjectHits(annotations)[queryHits(annotations) == i]
    unique(gr_reduced$repeat_class[hits])
  })
  
  # Now you can turn off the parallel plan
  plan(sequential)
  
  print(Sys.time())
  print("Repeat class assignment completed!")
  pb$tick()
  
  
  # Paste together annotations for each partitioned interval
  disjoined$annotation <- sapply(disjoined_annotations, function(annots) {
    paste(sort(annots), collapse = "-")
  })
  
  # Get seqinfo from the genome fasta
  seqInfo <- 
    as.data.frame(seqinfo(genome)) %>% 
    rownames_to_column("seqnames") %>% 
    select(c("seqnames", "seqlengths"))
  
  if(species == "Dmel"){
    seqInfo <- 
      seqInfo %>% 
      mutate(seqnames = str_remove(seqnames, "Drosophila melanogaster ")) %>% 
      mutate(seqnames = str_remove(seqnames, "sequence")) %>% 
      mutate(seqnames = str_replace(seqnames, "chromosome ", "chr")) %>% 
      arrange(seqnames)
  } else if(species == "Drho"){
    contig_names <- 
      unlist(regmatches(seqInfo$seqnames, 
                        gregexpr("contig_\\d+", seqInfo$seqnames, perl = TRUE)))
    NW_names <- 
      unlist(regmatches(seqInfo$seqnames, 
                        gregexpr("^\\w+[.]\\d+", seqInfo$seqnames, perl = TRUE)))
    
    seqInfo <- 
      seqInfo %>% 
      mutate(seqnames = paste(contig_names, NW_names, sep = " ")) %>% 
      arrange(seqnames)
  } else{
    seqInfo <- 
      seqInfo %>% 
      arrange(seqnames)
  }
  
  # Check seqnames of seqInfo and disjoined
  if(sum(levels(seqnames(disjoined)) %in% seqInfo$seqnames) == 0){
    matches <- 
      sapply(levels(seqnames(disjoined)), function(seqname){
        grep(seqname, seqInfo$seqnames)
      })
    num_matches <- sapply(matches, length)
    
    if(any(num_matches > 1)){
      msg <- sprintf("The following seqnames have multiple hits: %s", names(num_matches > 1))
      print(msg)
      stop("Multiple hits when mapping disjoined seqnames to seqInfo seqnames!")
    }
    else{
      seqlevels(disjoined) <- seqInfo$seqnames[matches]
    }
  } else if(sum(levels(seqnames(disjoined)) %in% seqInfo$seqnames) < length(levels(seqnames(disjoined)))){
    msg <- sprintf("The following seqnames have no match in seqInfo seqnames: %s", 
                   levels(seqnames(disjoined))[!levels(seqnames(disjoined)) %in% seqInfo$seqnames])
    warning(msg)
  }
  
  # Summarize nonoverlapping width by repeatClass
  repeatClassWidth <- 
    data.frame(seqnames = seqnames(disjoined), 
               repeatClass = disjoined$annotation,
               repeatWidth = width(disjoined)) %>% 
    group_by(seqnames, repeatClass) %>% 
    summarize(nonOverlapWidth = sum(repeatWidth, na.rm = TRUE)) %>% 
    spread(key = "repeatClass", value = "nonOverlapWidth") %>% 
    full_join(seqInfo, by = "seqnames") %>% 
    column_to_rownames("seqnames") %>% 
    mutate(nonrepeat = seqlengths - rowSums(select(., -seqlengths), na.rm = TRUE)) %>% 
    select(seqlengths, nonrepeat, everything()) %>% 
    rownames_to_column("seqnames")
  
  # Summarize nonoverlapping percentage by repeatClass
  repeatClassPerc <- 
    repeatClassWidth %>% 
    column_to_rownames("seqnames") %>% 
    mutate(across(-seqlengths, ~ . / seqlengths * 100)) %>% 
    rownames_to_column("seqnames")
  
  print("Writing output tables to CSV!")
  pb$tick()
  
  # Write detailed repeat content by scaffold to CSV
  write.csv(repeatClassWidth, file.path(outdata_dir, "scaffold_distribution_of_repeat_width.csv"), row.names = FALSE)
  write.csv(repeatClassPerc, file.path(outdata_dir, "scaffold_distribution_of_repeat_percentage.csv"), row.names = FALSE)
  
  # target repeat categories
  target_cats <- 
    c("DNA", "RC", "LINE", "SINE", "LTR", "tandem_repeat", "Unknown",
      "Simple", "rRNA", "tRNA", "snRNA", "nonrepeat", "Satellite")
  target_cats <-
    unique(disjoined$annotation)[unique(disjoined$annotation) %in% target_cats]
  
  # Create a stacked bar plot
  output <- 
    repeatClassPerc %>% 
    arrange(seqnames) %>% 
    column_to_rownames("seqnames") %>% 
    mutate(multiclass = 100 - rowSums(select(., nonrepeat, all_of(target_cats)),
                                      na.rm = TRUE)) %>% 
    mutate(multiclass = ifelse(multiclass < 0, 0, multiclass)) %>% 
    select(seqlengths, nonrepeat, multiclass, everything()) %>% 
    rownames_to_column("seqnames") %>% 
    gather(key = "repeatClass", value = "nonOverlapPerc", -seqnames, -seqlengths) %>% 
    filter(repeatClass %in% c("multiclass", "nonrepeat", target_cats)) %>% 
    # mutate(repeatClass = factor(repeatClass, 
    #                             levels = c("DNA", "Simple", "tandem_repeat", "RC", "LINE", "LTR", 
    #                                        "rRNA", "tRNA", "multiclass", "nonrepeat", "Unknown", 
    #                                        "SINE", "snRNA", "Satellite")),
    #        nonOverlapPerc = replace_na(nonOverlapPerc, 0))  
    mutate(repeatClass = factor(repeatClass),
           nonOverlapPerc = replace_na(nonOverlapPerc, 0))  
  
  # customize stacked bar plot by using geom_rect()
  output <- 
    output %>% 
    split(.$seqnames)  %>%
    lapply(function(df){
      ans <- 
        df %>% 
        mutate(ymin = cumsum(nonOverlapPerc) - nonOverlapPerc,
               ymax = cumsum(nonOverlapPerc))
      ans
    }) %>% 
    bind_rows()
  
  # update seqInfo with increasing scaffold length
  output <- 
    seqInfo %>% 
    arrange(seqlengths) %>% 
    mutate(xmin = sqrt(cumsum(seqlengths) - seqlengths),
           xmax = sqrt(cumsum(seqlengths))) %>% 
    select(c("seqnames", "xmin", "xmax")) %>% 
    full_join(output, ., by = "seqnames")
  
  # prepare summary tables
  output_tab <- 
    output %>% 
    group_by(repeatClass) %>% 
    summarize(width = sum(seqlengths * nonOverlapPerc / 100, na.rm = TRUE)) %>% 
    mutate(perc = 100 * width / sum(width, na.rm = TRUE))
    
  
  # Write detailed repeat content by scaffold to CSV
  write.csv(output_tab, file.path(outdata_dir, "repeat_content_by_class_summary.csv"), row.names = FALSE)
  
  # set text for x axis (suppress crowded text)
  xlab_text <- 
    output %>% 
    filter(xmax - xmin >= 100) %>% 
    mutate(x = (xmax + xmin)/2) %>% 
    select(c("seqnames", "x")) %>% 
    distinct()
  
  # Define color palettes
  orange_palette <- brewer.pal(3, "Oranges")
  green_palette <- brewer.pal(2, "Greens")[1:2]
  blue_palette <- brewer.pal(3, "Blues")
  purple_palette <- brewer.pal(3, "Purples")
  black_palette <- c("black")
  dark_gray_palette <- c("gray40")  # Dark gray
  light_gray_palette <- c("gray80")  # Light gray
  
  # Define the groups and assign the colors
  groups <- c("tandem_repeat", "Simple", "Satellite", 
              "DNA", "RC", 
              "LINE", "SINE", "LTR",
              "rRNA", "tRNA", "snRNA", 
              "Unknown", 
              "multiclass", 
              "nonrepeat")
  
  fills <- c(orange_palette, 
             green_palette, 
             blue_palette, 
             purple_palette, 
             black_palette, 
             dark_gray_palette, 
             light_gray_palette)
  
  # Create a named vector that matches groups to colors
  fill_mapping <- setNames(fills, groups)
  
  plot <- 
    output %>% 
    mutate(repeatClass = factor(repeatClass, levels = groups)) %>% 
    ggplot(aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = repeatClass)) +
    geom_rect(color = "black", lwd = 0.1, alpha = 0.75) +
    labs(x = ">1Mb scaffold shown", y = "Proportion (%)", fill = "") +
    scale_fill_manual(values = fill_mapping) + 
    scale_x_continuous(breaks = xlab_text$x, labels = xlab_text$seqnames) +
    guides(fill = guide_legend(override.aes = list(alpha = 0.75),
                               nrow = 2, byrow = TRUE),) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
  
  print("Saving output figures to PDF/PNG!")
  pb$tick()
  
  # Save plot
  ggsave(filename = file.path(outresult_dir, "scaffold_distribution_of_repeats.pdf"), 
         plot = plot, height = 6, width = 8)
  ggsave(filename = file.path(outresult_dir, "scaffold_distribution_of_repeats.png"), 
         plot = plot, height = 6, width = 8)
  
  pb$tick()
  return(disjoined)
}

# debugonce(processGenome)
# genomeFile <- "E:/Sequencing/Novogene/C202SC19030542/2020-11-18/References/prolongata_renamed_assembly.fasta"
# species <- "Dpro"
# bedFile1 <- "data/final/Dpro/RepeatModeler/renamed/families-classified.bed"
# bedFile2 <- "data/final/Dpro/TandemRepeatFinder/renamed/prolongata_renamed_assembly_TRF.bed"
# disjoined <- processGenome(genomeFile, "Dpro", bedFile1, bedFile2)

disjoined <- processGenome(genomeFile, species, bedFile1, bedFile2)
