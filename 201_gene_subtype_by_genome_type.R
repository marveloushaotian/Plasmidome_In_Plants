#!/usr/bin/env Rscript

# =============================================================================
# Top-N Gene Subtype Plots by Host and Contig_Type3
# Description: Plot top N subtypes as grouped bars for
#              Defense_Subtype, AntiDS_Type, and AMR_Type.
#              For each AMR subtype, show four hosts together.
#              Generate separate plots for each Contig_Type3.
# Usage: Rscript 201_gene_subtype_by_genome_type.R -i <input.csv> -o <output_dir> [-n <top_n>]
#
# Arguments:
#   -i: Input CSV file path (required)
#   -o: Output directory for results (required)
#   -n: Number of top subtypes to analyze for each type (default: 15)
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
parser <- ArgumentParser(description = "Generate statistics and plots for top gene subtypes")
parser$add_argument("-i", "--input", required = TRUE, 
                    help = "Input CSV file path (e.g., Contig_Sample_Mapping_Final.csv)")
parser$add_argument("-o", "--output", required = TRUE,
                    help = "Output directory for results")
parser$add_argument("-n", "--top-n", type = "integer", default = 15,
                    help = "Number of top subtypes to analyze for each type (default: 15)")

args <- parser$parse_args()

# Read CSV file
cat(sprintf("Reading data from: %s\n", args$input))
data <- read.csv(args$input, stringsAsFactors = FALSE, check.names = FALSE)

# Create output directory if it doesn't exist
if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", args$output))
}

# Validate required columns
required_cols <- c("Host", "Sample_ID", "Contig_Type3", "chromosome_length", "plasmid_length", "virus_length")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

# Build total kb per Host and Contig_Type3 using unique Sample_ID
length_by_host_ct3 <- data %>%
  select(Host, Sample_ID, chromosome_length, plasmid_length, virus_length) %>%
  mutate(Host = trimws(Host)) %>%
  filter(!is.na(Host), Host != "", !is.na(Sample_ID), Sample_ID != "") %>%
  distinct() %>%
  tidyr::pivot_longer(
    cols = c(chromosome_length, plasmid_length, virus_length),
    names_to = "length_type",
    values_to = "total_length_bp"
  ) %>%
  mutate(
    Contig_Type3 = dplyr::recode(
      length_type,
      "chromosome_length" = "Chromosome",
      "plasmid_length" = "Plasmid",
      "virus_length" = "Virus"
    )
  ) %>%
  group_by(Host, Contig_Type3) %>%
  summarise(total_length_kb = sum(total_length_bp, na.rm = TRUE) / 1000, .groups = "drop")

# Defense-related subtype names to exclude
defense_exclude_list <- c("VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA")

# Process one subtype column and generate plots
process_one_type <- function(df, col_name, output_prefix, top_n, defense_filter = FALSE) {
  if (!(col_name %in% colnames(df))) {
    cat(sprintf("Skipping %s: column not found.\n", col_name))
    return(invisible(NULL))
  }

  subtype_counts <- df %>%
    select(Host, Contig_Type3, all_of(col_name)) %>%
    mutate(
      Host = trimws(Host),
      Contig_Type3 = dplyr::if_else(trimws(Contig_Type3) %in% c("Virus", "Provirus"), "Virus", trimws(Contig_Type3)),
      subtype = strsplit(.[[col_name]], ",")
    ) %>%
    unnest(subtype) %>%
    mutate(subtype = trimws(subtype)) %>%
    filter(
      !is.na(Host), Host != "",
      !is.na(Contig_Type3), Contig_Type3 != "",
      !is.na(subtype), subtype != ""
    )

  if (defense_filter) {
    subtype_counts <- subtype_counts %>%
      filter(!(subtype %in% defense_exclude_list), tolower(subtype) != "other")
  }

  subtype_counts <- subtype_counts %>%
    count(Host, Contig_Type3, subtype, name = "count")

  host_order <- subtype_counts %>%
    group_by(Host) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    arrange(desc(total_count)) %>%
    slice_head(n = 4) %>%
    pull(Host)

  subtype_counts <- subtype_counts %>%
    filter(Host %in% host_order)

  top_subtypes <- subtype_counts %>%
    group_by(subtype) %>%
    summarise(total_count = sum(count), .groups = "drop") %>%
    arrange(desc(total_count)) %>%
    slice_head(n = top_n) %>%
    pull(subtype)

  subtype_counts <- subtype_counts %>%
    filter(subtype %in% top_subtypes)

  contig_type3_order <- sort(unique(subtype_counts$Contig_Type3))
  all_combos <- expand.grid(
    Host = host_order,
    Contig_Type3 = contig_type3_order,
    subtype = top_subtypes,
    stringsAsFactors = FALSE
  )

  counts_top_full <- all_combos %>%
    left_join(subtype_counts, by = c("Host", "Contig_Type3", "subtype")) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    left_join(length_by_host_ct3, by = c("Host", "Contig_Type3")) %>%
    mutate(
      count_perkb = ifelse(
        !is.na(total_length_kb) & total_length_kb > 0,
        count / total_length_kb,
        NA_real_
      )
    )

  subtype_levels <- counts_top_full %>%
    group_by(subtype) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(subtype)

  counts_top_full <- counts_top_full %>%
    mutate(
      subtype = factor(subtype, levels = subtype_levels),
      Host = factor(Host, levels = host_order),
      Contig_Type3 = factor(Contig_Type3, levels = contig_type3_order)
    )

  csv_file <- file.path(args$output, paste0(output_prefix, "_Top", top_n, "_CountPerkb_ByHost_ContigType3.csv"))
  write_csv(counts_top_full, csv_file)
  cat(sprintf("Table saved to: %s\n", csv_file))
  cat(sprintf("Selected hosts (top 4): %s\n", paste(host_order, collapse = ", ")))
  cat(sprintf("Selected top %d subtypes for %s.\n", length(top_subtypes), col_name))

  for (ct3 in contig_type3_order) {
    plot_df <- counts_top_full %>%
      filter(Contig_Type3 == ct3)

    p <- ggplot(plot_df, aes(x = subtype, y = count_perkb, fill = Host)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.75) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)
      ) +
      labs(
        title = paste0("Top ", length(top_subtypes), " ", col_name, " - ", ct3),
        x = col_name,
        y = "Count per kb",
        fill = "Host"
      )

    ct3_file_tag <- gsub("[^A-Za-z0-9]+", "_", ct3)
    pdf_file <- file.path(
      args$output,
      paste0(output_prefix, "_Top", top_n, "_CountPerkb_ByHost_", ct3_file_tag, ".pdf")
    )
    ggsave(pdf_file, plot = p, width = 16, height = 7)
    cat(sprintf("Plot saved to: %s\n", pdf_file))
  }
}

cat("\n========================================\n")
cat("Processing Defense_Subtype...\n")
cat("========================================\n")
process_one_type(data, "Defense_Subtype", "Defense_Subtype", args$top_n, defense_filter = TRUE)

cat("\n========================================\n")
cat("Processing AntiDS_Type...\n")
cat("========================================\n")
process_one_type(data, "AntiDS_Type", "AntiDS_Type", args$top_n, defense_filter = TRUE)

cat("\n========================================\n")
cat("Processing AMR_Type...\n")
cat("========================================\n")
process_one_type(data, "AMR_Type", "AMR_Type", args$top_n, defense_filter = FALSE)

cat("\n========================================\n")
cat("All processing completed!\n")
cat("========================================\n")
cat(sprintf("\nOutput directory: %s\n", args$output))
cat(sprintf("Output files:\n"))
cat(sprintf("  - Defense_Subtype_Top%d_CountPerkb_ByHost_ContigType3.csv\n", args$top_n))
cat(sprintf("  - Defense_Subtype_Top%d_CountPerkb_ByHost_<Contig_Type3>.pdf\n", args$top_n))
cat(sprintf("  - AntiDS_Type_Top%d_CountPerkb_ByHost_ContigType3.csv\n", args$top_n))
cat(sprintf("  - AntiDS_Type_Top%d_CountPerkb_ByHost_<Contig_Type3>.pdf\n", args$top_n))
cat(sprintf("  - AMR_Type_Top%d_CountPerkb_ByHost_ContigType3.csv\n", args$top_n))
cat(sprintf("  - AMR_Type_Top%d_CountPerkb_ByHost_<Contig_Type3>.pdf\n", args$top_n))
cat("\n")

