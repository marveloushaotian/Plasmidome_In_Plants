#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Step 1: Set input and output paths
input_file <- "Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final.csv"
output_dir <- "Result/NCBI_4395_Batch/03_Alpha_Diversity"

# Step 2: Create output directory if not exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Step 3: Read input data
cat("Reading input data...\n")
data <- fread(input_file, header = TRUE, stringsAsFactors = FALSE)
cat(sprintf("Loaded %d rows and %d columns\n", nrow(data), ncol(data)))

# Step 4: Function to split comma-separated values and count unique items
count_unique_types <- function(values) {
  # Remove NA, empty strings, and whitespace
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) return(0)
  
  # Split by comma and trim whitespace
  split_values <- unlist(strsplit(values, ","))
  split_values <- trimws(split_values)
  
  # Remove empty strings after splitting
  split_values <- split_values[split_values != ""]
  
  # Count unique non-empty values
  return(length(unique(split_values)))
}

# Step 5: Calculate richness for Defense_Type
cat("Calculating Defense_Type richness...\n")
defense_richness <- data %>%
  filter(!is.na(Defense_Type) & Defense_Type != "") %>%
  group_by(Sample_ID, Host, Contig_Type) %>%
  summarise(
    Defense_Richness = count_unique_types(Defense_Type),
    .groups = "drop"
  )

# Step 6: Calculate richness for AntiDS_Type
cat("Calculating AntiDS_Type richness...\n")
antids_richness <- data %>%
  filter(!is.na(AntiDS_Type) & AntiDS_Type != "") %>%
  group_by(Sample_ID, Host, Contig_Type) %>%
  summarise(
    AntiDS_Richness = count_unique_types(AntiDS_Type),
    .groups = "drop"
  )

# Step 7: Calculate richness for AMR_Type
cat("Calculating AMR_Type richness...\n")
amr_richness <- data %>%
  filter(!is.na(AMR_Type) & AMR_Type != "") %>%
  group_by(Sample_ID, Host, Contig_Type) %>%
  summarise(
    AMR_Richness = count_unique_types(AMR_Type),
    .groups = "drop"
  )

# Step 8: Merge all richness data
cat("Merging richness data...\n")
all_samples <- data %>%
  select(Sample_ID, Host, Contig_Type) %>%
  distinct()

richness_data <- all_samples %>%
  left_join(defense_richness, by = c("Sample_ID", "Host", "Contig_Type")) %>%
  left_join(antids_richness, by = c("Sample_ID", "Host", "Contig_Type")) %>%
  left_join(amr_richness, by = c("Sample_ID", "Host", "Contig_Type")) %>%
  mutate(
    Defense_Richness = ifelse(is.na(Defense_Richness), 0, Defense_Richness),
    AntiDS_Richness = ifelse(is.na(AntiDS_Richness), 0, AntiDS_Richness),
    AMR_Richness = ifelse(is.na(AMR_Richness), 0, AMR_Richness)
  )

# Step 9: Save richness data
output_file <- file.path(output_dir, "defense_richness_by_sample.csv")
cat(sprintf("Saving richness data to %s...\n", output_file))
fwrite(richness_data, output_file)

# Step 10: Prepare data for plotting (long format)
richness_long <- richness_data %>%
  pivot_longer(
    cols = c(Defense_Richness, AntiDS_Richness, AMR_Richness),
    names_to = "Type",
    values_to = "Richness"
  ) %>%
  mutate(
    Type = gsub("_Richness", "", Type)
  )

# Step 11: Plot Defense_Type richness by Host and Contig_Type
cat("Generating Defense_Type richness boxplot...\n")
p1 <- ggplot(richness_data, aes(x = Contig_Type, y = Defense_Richness, fill = Host)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Host, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Defense System Richness per Sample",
    subtitle = "Grouped by Host and Contig Type",
    x = "Contig Type",
    y = "Defense System Richness"
  )

pdf_file <- file.path(output_dir, "Defense_Richness_by_Host_ContigType.pdf")
ggsave(pdf_file, p1, width = 12, height = 10)
cat(sprintf("Saved plot to %s\n", pdf_file))

# Step 12: Plot AntiDS_Type richness by Host and Contig_Type
cat("Generating AntiDS_Type richness boxplot...\n")
p2 <- ggplot(richness_data, aes(x = Contig_Type, y = AntiDS_Richness, fill = Host)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Host, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Anti-Defense System Richness per Sample",
    subtitle = "Grouped by Host and Contig Type",
    x = "Contig Type",
    y = "Anti-Defense System Richness"
  )

pdf_file <- file.path(output_dir, "AntiDS_Richness_by_Host_ContigType.pdf")
ggsave(pdf_file, p2, width = 12, height = 10)
cat(sprintf("Saved plot to %s\n", pdf_file))

# Step 13: Plot AMR_Type richness by Host and Contig_Type
cat("Generating AMR_Type richness boxplot...\n")
p3 <- ggplot(richness_data, aes(x = Contig_Type, y = AMR_Richness, fill = Host)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~Host, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "AMR Type Richness per Sample",
    subtitle = "Grouped by Host and Contig Type",
    x = "Contig Type",
    y = "AMR Type Richness"
  )

pdf_file <- file.path(output_dir, "AMR_Richness_by_Host_ContigType.pdf")
ggsave(pdf_file, p3, width = 12, height = 10)
cat(sprintf("Saved plot to %s\n", pdf_file))

# Step 14: Create combined plot for all three types
cat("Generating combined richness boxplot...\n")
p4 <- ggplot(richness_long, aes(x = Contig_Type, y = Richness, fill = Type)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(Type ~ Host, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Defense, Anti-Defense, and AMR System Richness per Sample",
    subtitle = "Grouped by Host and Contig Type",
    x = "Contig Type",
    y = "Richness",
    fill = "System Type"
  )

pdf_file <- file.path(output_dir, "Combined_Richness_by_Host_ContigType.pdf")
ggsave(pdf_file, p4, width = 14, height = 12)
cat(sprintf("Saved plot to %s\n", pdf_file))

# Step 15: Generate summary statistics
cat("Generating summary statistics...\n")
summary_stats <- richness_data %>%
  group_by(Host, Contig_Type) %>%
  summarise(
    n_samples = n(),
    Defense_mean = mean(Defense_Richness, na.rm = TRUE),
    Defense_median = median(Defense_Richness, na.rm = TRUE),
    Defense_sd = sd(Defense_Richness, na.rm = TRUE),
    AntiDS_mean = mean(AntiDS_Richness, na.rm = TRUE),
    AntiDS_median = median(AntiDS_Richness, na.rm = TRUE),
    AntiDS_sd = sd(AntiDS_Richness, na.rm = TRUE),
    AMR_mean = mean(AMR_Richness, na.rm = TRUE),
    AMR_median = median(AMR_Richness, na.rm = TRUE),
    AMR_sd = sd(AMR_Richness, na.rm = TRUE),
    .groups = "drop"
  )

summary_file <- file.path(output_dir, "richness_summary_statistics.csv")
fwrite(summary_stats, summary_file)
cat(sprintf("Saved summary statistics to %s\n", summary_file))

cat("\nAnalysis completed successfully!\n")
cat(sprintf("Results saved in: %s\n", output_dir))

