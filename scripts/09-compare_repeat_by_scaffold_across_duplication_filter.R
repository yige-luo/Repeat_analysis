# Load libraries
if (!requireNamespace("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
}
if (!requireNamespace("patchwork", quietly = TRUE)){
  install.packages("patchwork")
}
library(tidyverse)
library(patchwork)

# species under analysis
species <- c("Dpro")
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

if(!dir.exists(outresult_dir)){
  dir.create(outresult_dir, recursive = TRUE, showWarnings = FALSE)
}


# Load data
repeat_df <- read_csv(file.path(input_dir, "scaffold_distribution_of_repeat_percentage.csv"))
annotation_df <- read_csv("data/final//Dpro/scaffold_duplication_annotation.csv")

# combine tables
merged_df <- 
  repeat_df %>% 
  select("seqnames", "seqlengths", "nonrepeat") %>% 
  mutate(`repeat` = 100 - nonrepeat) %>% 
  left_join(annotation_df, ., by = "seqnames")

plot1 <- 
  merged_df %>% 
  ggplot(aes(`repeat`, `seqlengths`, color = `filter`)) +
  geom_point() +
  labs(x = "", color = "") + 
  scale_y_log10(name = "Scaffold length (bp)") +
  facet_wrap(~filter, nrow = 1) +
  theme_bw() + 
  theme(strip.text = element_text(face = "bold", size = 12), 
        legend.position = "none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, -0.5, 0, "in"))
# plot1

plot2 <- 
  merged_df %>% 
  ggplot(aes(`repeat`, fill = `filter`)) +
  geom_histogram(aes(y = after_stat(density)), color = "black", 
                 binwidth = 10) +
  labs(x = "Repeat content (%)") +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  scale_y_continuous(breaks = c(0, 0.02, 0.04)) +
  facet_wrap(~filter, nrow = 1) +
  theme_minimal() + 
  theme(legend.position = "none", panel.grid = element_blank(),
        strip.text = element_blank(), strip.background = element_blank(), 
        axis.ticks.x = element_line(),
        # axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(-0.5, 0, 0, 0, "in"))
plot2

# Combine the two plots
combined_plot <- 
  plot1 / plot2 +
  plot_layout(heights = c(4, 1))

# # Print the combined plot
# print(combined_plot)
# Save plot
ggsave(filename = file.path(outresult_dir, "scaffold_repeat_profile_by_duplication.pdf"), 
       plot = combined_plot, height = 5, width = 7.5)
ggsave(filename = file.path(outresult_dir, "scaffold_repeat_profile_by_duplication.png"), 
       plot = combined_plot, height = 5, width = 7.5)
