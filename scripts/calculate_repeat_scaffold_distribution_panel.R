# Load libraries
library(tidyverse)
library(RColorBrewer)
library(grid)
library(gridExtra)

# Load data
args <- commandArgs(trailingOnly = TRUE)
repeatFile <- args[1]
annoFile <- args[2]
seqinfoFile <- args[3]
species <- args[4]
repeat_df <- read_csv(repeatFile)
annotation_df <- read_csv(annoFile)
seqinfo_df <- read_csv(seqinfoFile)

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

# species under analysis
# species <- "Dpro"
# Prepare for output
output_dir <- paste0("results/", species, "/", "Table_Plot_Prototype", "/", Sys.Date(), "/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
# input_date <- "2023-06-29a"
# input_dir <- paste0("results/", species, "/", "Table_Plot_Prototype", "/", input_date, "/")
# repeat_df <- read_csv(file.path(input_dir, "scaffold_distribution_of_repeat_percentage.csv"))
# annotation_df <- read_csv("data/scaffold_duplication_annotation.csv")
# seqinfo_df <- read_csv("data/Dpro_seqInfo.csv")

# Keep only the columns we want
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
         nonOverlapPerc = replace_na(nonOverlapPerc, 0))  

# update annotation_df
annotation_df <- 
  seqinfo_df %>% 
  right_join(annotation_df, by = "seqnames") %>% 
  mutate(seq_length_range = cut_number(seqlengths, n = ceiling(n()/36)))

# Merge data frames
merged_df <- full_join(select(repeat_df, -c("seqlengths")), 
                       annotation_df, by = "seqnames")

# Generate the base plot
p <- ggplot(data = merged_df, aes(x = seqnames, y = nonOverlapPerc, fill = repeatClass)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = fill_mapping) +
  facet_wrap(~seq_length_range, scales = "free_x") +
  labs(x = NULL, y = "Percentage", fill = NULL) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.75),
                           nrow = 2, byrow = TRUE),) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6, face = "bold"))
# p

# Extract the plot grob
g <- ggplotGrob(p)

# Get the names of the bottom axis grobs
axis_b_names <- grep("^axis-b", g$layout$name, value = TRUE)

# Define color mapping for the "filter" column
filter_colors <- c("remove" = "red", "keep" = "blue", "as_is" = "black")
filter_font <- c("remove" = 2L, "keep" = 2L, "as_is" = 1L)

# Iterate through each bottom axis grob
for (axis_b_name in axis_b_names) {
  # Extract the bottom axis grob
  xaxis_grob <- g$grobs[[which(g$layout$name == axis_b_name)]] 
  
  # Extract the tick labels from the bottom axis grob
  tick_labels <- xaxis_grob$children$axis$grobs[[2]]$children[[1]]$label
  
  # Create a data frame of the x-axis tick labels and their corresponding "filter" values
  tick_label_df <- 
    annotation_df %>% 
    filter(seqnames %in% tick_labels) %>% 
    arrange(seqnames)
  
  
  # Map the "filter" values to colors
  tick_label_df$color <- filter_colors[tick_label_df$filter]
  tick_label_df$font <- filter_font[tick_label_df$filter]
  
  # Modify the bottom axis grob to add colored tick labels
  xaxis_grob$children$axis$grobs[[2]]$children[[1]]$gp$col <- tick_label_df$color
  xaxis_grob$children$axis$grobs[[2]]$children[[1]]$gp$font <- tick_label_df$font
  
  # Replace the original bottom axis grob with the modified bottom axis grob
  g$grobs[[which(g$layout$name == axis_b_name)]] <- xaxis_grob
}

# Open PDF device
pdf(paste0(output_dir, "scaffold_distribution_of_repeats_panels.pdf"), width = 11, height = 6)

# Draw the plot
grid.draw(g)

# Close the PDF device
dev.off()

# Open PDF device
png(paste0(output_dir, "scaffold_distribution_of_repeats_panels.png"), 
    width = 11, height = 6, units = "in", res = 300)

# Draw the plot
grid.draw(g)

# Close the PDF device
dev.off()


