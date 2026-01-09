#!/usr/bin/env Rscript

# =============================================================================
# Statistics and Plots for Gene Subtypes
# Description: Generate statistics and bar plots for top N subtypes of 
#              Defense_Subtype, AntiDS_Type, and AMR_Type
#              Counts subtypes by Host and Contig_Type, filters defense-related
#              subtypes, and creates visualization
# Usage: Rscript 201_statistics_plot.R -i <input.csv> -o <output_dir> [-n <top_n>]
#
# Arguments:
#   -i: Input CSV file path (required)
#   -o: Output directory for results (required)
#   -n: Number of top subtypes to analyze (default: 100)
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Generate statistics and plots for gene subtypes")
parser$add_argument("-i", "--input", required = TRUE, 
                    help = "Input CSV file path (e.g., Contig_Sample_Mapping_Final.csv)")
parser$add_argument("-o", "--output", required = TRUE,
                    help = "Output directory for results")
parser$add_argument("-n", "--top-n", type = "integer", default = 100,
                    help = "Number of top subtypes to analyze (default: 100)")

args <- parser$parse_args()

# Read CSV file
cat(sprintf("Reading data from: %s\n", args$input))
data <- read.csv(args$input, stringsAsFactors = FALSE, check.names = FALSE)

# Create output directory if it doesn't exist
if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", args$output))
}

# Merge Contig_Type (replace Provirus with Virus)
chl <- dplyr::mutate(data,
                     Contig_Type = dplyr::if_else(Contig_Type %in% c('Virus', 'Provirus'), 'Virus', Contig_Type))

# Defense system subtypes to exclude (same as taxonomy_heatmap.R)
defense_exclude_list <- c("VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA")

# Define function for splitting and counting subtypes, with optional defense filter
# 额外: 对Defense_Subtype和AntiDS_Type需要过滤掉"Other"
count_subtypes <- function(df, col_name, defense_filter = FALSE) {
  res <- df %>%
    select(Host, Contig_Type, all_of(col_name)) %>%
    mutate(subtype = strsplit(.[[col_name]], ',')) %>%
    unnest(subtype) %>%
    mutate(subtype = trimws(subtype)) %>%
    filter(!is.na(subtype), subtype != "") 
  
  # 筛选规则：
  # - 防御类（defense_filter==TRUE）：去除defense_exclude_list和"Other"
  # - 非防御类：不做特殊过滤
  if (defense_filter) {
    res <- res %>% filter(!(subtype %in% defense_exclude_list), tolower(subtype) != "other")
  }
  res %>%
    count(Host, Contig_Type, subtype, sort = TRUE, name = 'count')
}

# Function to process and plot top N subtypes, all bars in gray, with optional defense filtering
process_and_plot <- function(df, col_name, output_prefix, output_dir, top_n, defense_filter = FALSE) {
  # Count subtypes
  subtype_counts <- count_subtypes(df, col_name, defense_filter)
  
  # Get top N subtypes based on total count
  top_subtypes <- subtype_counts %>%
    group_by(subtype) %>%
    summarise(total_count = sum(count), .groups = 'drop') %>%
    arrange(desc(total_count)) %>%
    slice_head(n = top_n) %>%
    pull(subtype)
  
  cat(sprintf("  Found %d unique subtypes, selecting top %d\n", 
              length(unique(subtype_counts$subtype)), length(top_subtypes)))
  
  # Keep only top N subtypes
  counts_topN <- subtype_counts %>%
    filter(subtype %in% top_subtypes)
  
  # Fill in missing combinations with zeros
  all_combos <- expand.grid(
    Host = unique(counts_topN$Host),
    Contig_Type = unique(counts_topN$Contig_Type),
    subtype = top_subtypes,
    stringsAsFactors = FALSE
  )
  
  counts_topN_full <- all_combos %>%
    left_join(counts_topN, by = c("Host", "Contig_Type", "subtype")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  # Set subtype as factor with descending order
  subtype_levels <- counts_topN %>%
    group_by(subtype) %>%
    summarise(total = sum(count), .groups = 'drop') %>%
    arrange(desc(total)) %>%
    pull(subtype)
  
  counts_topN_full$subtype <- factor(counts_topN_full$subtype, levels = subtype_levels)
  
  # Save table
  csv_file <- file.path(output_dir, paste0(output_prefix, '_Top', top_n, '_Counts_ByGroup.csv'))
  write_csv(counts_topN_full, csv_file)
  cat(sprintf("  Table saved to: %s\n", csv_file))
  
  # Plot (all gray)
  p <- ggplot(counts_topN_full, aes(x = subtype, y = count)) +
    geom_bar(stat = 'identity', position = 'dodge', fill = "gray60") +
    facet_wrap(~Host + Contig_Type, scales = "fixed", ncol = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    labs(title = paste0('Top ', top_n, ' ', col_name, ' Counts by Host and Contig Type'),
         x = col_name, y = 'Count')
  
  pdf_file <- file.path(output_dir, paste0(output_prefix, '_Top', top_n, '_Counts_ByGroup.pdf'))
  ggsave(
    pdf_file,
    plot = p,
    width = 20, height = 30,
    limitsize = FALSE
  )
  cat(sprintf("  Plot saved to: %s\n", pdf_file))
  
  return(counts_topN_full)
}

# Process Defense_Subtype (with defense filter, also remove "Other")
cat("\n========================================\n")
cat("Processing Defense_Subtype...\n")
cat("========================================\n")
defense_results <- process_and_plot(chl, 'Defense_Subtype', 'Defense_Subtype', 
                                     args$output, args$top_n, defense_filter = TRUE)

# Process AntiDS_Type (with defense filter, also remove "Other")
cat("\n========================================\n")
cat("Processing AntiDS_Type...\n")
cat("========================================\n")
antids_results <- process_and_plot(chl, 'AntiDS_Type', 'AntiDS_Type', 
                                     args$output, args$top_n, defense_filter = TRUE)

# Process AMR_Type (no filter)
cat("\n========================================\n")
cat("Processing AMR_Type...\n")
cat("========================================\n")
amr_results <- process_and_plot(chl, 'AMR_Type', 'AMR_Type', 
                                 args$output, args$top_n, defense_filter = FALSE)

cat("\n========================================\n")
cat("All processing completed!\n")
cat("========================================\n")
cat(sprintf("\nOutput directory: %s\n", args$output))
cat(sprintf("Output files:\n"))
cat(sprintf("  - Defense_Subtype_Top%d_Counts_ByGroup.csv/pdf\n", args$top_n))
cat(sprintf("  - AntiDS_Type_Top%d_Counts_ByGroup.csv/pdf\n", args$top_n))
cat(sprintf("  - AMR_Type_Top%d_Counts_ByGroup.csv/pdf\n", args$top_n))
cat("\n")

