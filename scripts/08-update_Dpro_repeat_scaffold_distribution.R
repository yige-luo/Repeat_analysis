# Load libraries
if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)){
  install.packages("RColorBrewer")
}
library(tidyverse)
library(RColorBrewer)

# Load data

# species under analysis
species <- "Dpro"
analysis <- "Table_Plot_Prototype"
partition <- ""
if(species == "Dpro"){
  # partition <- "removed"
  # partition <- "deduped"
  partition <- "renamed"
}
input_date <- "2023-07-20"
input_dir <- file.path("data/interim/", species, analysis, partition, input_date)

# Prepare for output
outresult_dir <- file.path("results", species, analysis, partition, Sys.Date())
outdata_dir <- file.path("data/interim", species, analysis, partition, Sys.Date())

if(!dir.exists(outresult_dir)){
  dir.create(outresult_dir, recursive = TRUE, showWarnings = FALSE)
}
if(!dir.exists(outdata_dir)){
  dir.create(outdata_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load data
repeat_df <- read_csv(file.path(input_dir, "scaffold_distribution_of_repeat_percentage.csv"))
annotation_df <- read_csv("data/final/Dpro/scaffold_duplication_annotation.csv")
seqinfo_df <- read_csv("data/final/Dpro/Dpro_seqInfo.csv")


# Define the groups and assign the colors
groups <- c("tandem_repeat", "Simple", "Satellite", 
            "DNA", "RC", 
            "LINE", "SINE", "LTR",
            "rRNA", "tRNA", "snRNA", 
            "Unknown", 
            "multiclass", 
            "nonrepeat")

# Keep only the columns we want, keep retained scaffolds
repeat_df <-    
  repeat_df %>% 
  arrange(seqnames) %>% 
  column_to_rownames("seqnames") %>% 
  mutate(multiclass = 100 - rowSums(select(., nonrepeat, one_of(groups)),
                                    na.rm = TRUE)) %>% 
  mutate(multiclass = ifelse(multiclass < 0, 0, multiclass)) %>% 
  select(seqlengths, nonrepeat, multiclass, everything()) %>% 
  rownames_to_column("seqnames") %>% 
  gather(key = "repeatClass", value = "nonOverlapPerc", -seqnames, -seqlengths) %>% 
  filter(repeatClass %in% c("multiclass", "nonrepeat", groups)) %>% 
  mutate(repeatClass = factor(repeatClass),
         nonOverlapPerc = replace_na(nonOverlapPerc, 0)) %>% 
  left_join(annotation_df, by = "seqnames") %>% 
  filter(filter != "remove")


# customize stacked bar plot by using geom_rect()
repeat_df <- 
  repeat_df %>% 
  split(.$seqnames)  %>%
  lapply(function(df){
    ans <- 
      df %>% 
      mutate(ymin = cumsum(nonOverlapPerc) - nonOverlapPerc,
             ymax = cumsum(nonOverlapPerc))
    ans
  }) %>% 
  bind_rows()

repeat_df <- 
  seqinfo_df %>% 
  arrange(seqlengths) %>% 
  mutate(xmin = sqrt(cumsum(seqlengths) - seqlengths),
         xmax = sqrt(cumsum(seqlengths))) %>% 
  select(c("seqnames", "xmin", "xmax")) %>% 
  right_join(., repeat_df, by = "seqnames")

# prepare summary tables
output_tab <- 
  repeat_df %>% 
  group_by(repeatClass) %>% 
  summarize(width = sum(seqlengths * nonOverlapPerc / 100, na.rm = TRUE)) %>% 
  mutate(perc = 100 * width / sum(width, na.rm = TRUE))

# Write detailed repeat content by scaffold to CSV
write.csv(output_tab, file.path(outdata_dir, "repeat_content_by_class_summary.csv"), row.names = FALSE)

# set text for x axis (suppress crowded text)
xlab_text <- 
  repeat_df %>% 
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
  repeat_df %>% 
  mutate(repeatClass = factor(repeatClass, levels = groups)) %>% 
  ggplot(aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax, fill = repeatClass)) +
  geom_rect(color = "black", lwd = 0.1, alpha = 0.75) +
  labs(x = "", y = "Proportion (%)", fill = "") +
  scale_fill_manual(values = fill_mapping) + 
  scale_x_continuous(breaks = xlab_text$x, labels = xlab_text$seqnames) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.75),
                             nrow = 2, byrow = TRUE),) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))

# Save plot
ggsave(filename = file.path(outresult_dir, "scaffold_distribution_of_repeats.pdf"), 
       plot = plot, height = 6, width = 8)
ggsave(filename = file.path(outresult_dir, "scaffold_distribution_of_repeats.png"), 
       plot = plot, height = 6, width = 8)
