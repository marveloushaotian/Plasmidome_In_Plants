#!/usr/bin/env Rscript

# =============================================================================
# Step 1: PCoA Calculation (Generic)
# Description: Calculate PCoA coordinates from gene type profiles
#              Supports: Defense, AMR, AntiDefense
#              Save coordinates and metadata for later plotting
# Usage: Rscript 213_step1_pcoa_calculate.R -i <input.csv> -t <type> -o <output_prefix> [-g <group_column>] [-c]
#
# Arguments:
#   -i: Input CSV file (expanded gene type file)
#   -t: Gene type: 'defense', 'amr', or 'antidefense'
#   -o: Output prefix (default: PCoA_Data)
#   -g: Optional grouping column (e.g., Host, Genus_CRBC_Updated) to split analysis by group
#   -c: Use count numbers instead of binary (present/absent) for PCoA analysis
# =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(ape)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Calculate PCoA coordinates from gene type profiles")
parser$add_argument("-i", "--input", required = TRUE, help = "Input CSV file path")
parser$add_argument("-t", "--type", required = TRUE, 
                     choices = c("defense", "amr", "antidefense"),
                     help = "Gene type: 'defense', 'amr', or 'antidefense'")
parser$add_argument("-o", "--output", default = "PCoA_Data", 
                    help = "Output file prefix (default: PCoA_Data)")
parser$add_argument("-g", "--group", default = NULL,
                    help = "Optional grouping column name to split analysis (e.g., Host, Genus_CRBC_Updated)")
parser$add_argument("-c", "--use-counts", action = "store_true",
                    help = "Use count numbers instead of binary (present/absent) for PCoA analysis")

args <- parser$parse_args()

# Configuration based on type
type_config <- list(
  defense = list(
    prefix = "^Defense_",
    exclude = c("Defense_VSPR", "Defense_dXTPase", "Defense_HEC-01", 
                "Defense_HEC-09", "Defense_PifA"),
    name = "Defense",
    matrix_suffix = "defense_matrix"
  ),
  amr = list(
    prefix = "^AMR_",
    exclude = character(0),
    name = "AMR",
    matrix_suffix = "amr_matrix"
  ),
  antidefense = list(
    prefix = "^AntiDS_",
    exclude = character(0),
    name = "AntiDefense",
    matrix_suffix = "antids_matrix"
  )
)

config <- type_config[[args$type]]

cat(sprintf("=== Step 1: PCoA Calculation for %s ===\n\n", config$name))
cat(sprintf("Reading data from: %s\n", args$input))
if (args$use_counts) {
  cat("Mode: Using count numbers (Bray-Curtis distance)\n")
} else {
  cat("Mode: Using binary presence/absence (Jaccard distance)\n")
}
cat("\n")
df <- read.csv(args$input, stringsAsFactors = FALSE, check.names = FALSE)

# Get input file directory for output files
input_dir <- dirname(normalizePath(args$input, mustWork = TRUE))

# Get numeric gene type columns
all_cols <- colnames(df)
gene_like <- all_cols[grepl(config$prefix, all_cols)]
existing_gene_cols <- gene_like[sapply(df[1:min(100, nrow(df)), gene_like], function(x) is.numeric(x))]

cat(sprintf("Numeric %s columns found (before filtering): %d\n", config$name, length(existing_gene_cols)))

# Apply exclusions
if (length(config$exclude) > 0) {
  existing_gene_cols <- existing_gene_cols[!existing_gene_cols %in% config$exclude]
  if (args$type == "defense") {
    existing_gene_cols <- existing_gene_cols[!grepl("other", existing_gene_cols, ignore.case = TRUE)]
  }
}

cat(sprintf("Numeric %s columns found (after filtering): %d\n", config$name, length(existing_gene_cols)))

# Filter for Plasmid and Chromosome only
df_filtered <- df %>%
  filter(Contig_Type2 %in% c("Plasmid", "Chromosome"))

cat(sprintf("Total rows after filtering: %d\n\n", nrow(df_filtered)))

# Check if grouping column exists
if (!is.null(args$group)) {
  if (!args$group %in% colnames(df_filtered)) {
    stop(sprintf("Error: Grouping column '%s' not found in input file. Available columns: %s\n", 
                 args$group, paste(head(colnames(df_filtered), 20), collapse = ", ")))
  }
  # Check for missing values in grouping column
  df_filtered <- df_filtered[!is.na(df_filtered[[args$group]]) & df_filtered[[args$group]] != "", ]
  groups <- unique(df_filtered[[args$group]])
  cat(sprintf("Grouping by '%s': Found %d groups\n", args$group, length(groups)))
  cat(sprintf("Groups: %s\n\n", paste(groups, collapse = ", ")))
} else {
  groups <- NULL
}

# Function to create gene matrix aggregated by Sample_ID
create_gene_matrix_by_sample <- function(data, contig_type, gene_cols) {
  data_sub <- data %>% filter(Contig_Type2 == contig_type)
  
  if (nrow(data_sub) == 0) return(NULL)
  
  # Aggregate by Sample_ID with all metadata
  agg_data <- data_sub %>%
    group_by(Sample_ID) %>%
    summarise(across(all_of(gene_cols), \(x) sum(x, na.rm = TRUE)),
              Host = first(Host),
              Class_CRBC = first(Class_CRBC),
              Order_CRBC = first(Order_CRBC),
              Family_CRBC = first(Family_CRBC),
              Genus_CRBC = first(Genus_CRBC),
              .groups = "drop")
  
  gene_matrix <- as.matrix(agg_data[, gene_cols])
  rownames(gene_matrix) <- agg_data$Sample_ID
  
  sample_info <- agg_data %>%
    select(Sample_ID, Host, Class_CRBC, Order_CRBC, Family_CRBC, Genus_CRBC)
  
  return(list(matrix = gene_matrix, sample_info = sample_info))
}

# Function to run PCoA and save results
run_pcoa_calculation <- function(data, contig_type, output_prefix, gene_cols, matrix_suffix, output_dir) {
  cat(sprintf("=== Processing %s ===\n", contig_type))
  
  result <- create_gene_matrix_by_sample(data, contig_type, gene_cols)
  if (is.null(result)) {
    cat(sprintf("No data for %s\n\n", contig_type))
    return(NULL)
  }
  
  gene_matrix <- result$matrix
  sample_info <- result$sample_info
  
  cat(sprintf("Number of samples: %d\n", nrow(gene_matrix)))
  
  # Remove samples with no genes
  row_sums <- rowSums(gene_matrix)
  gene_matrix_nonzero <- gene_matrix[row_sums > 0, , drop = FALSE]
  sample_info_nonzero <- sample_info %>% filter(Sample_ID %in% rownames(gene_matrix_nonzero))
  
  cat(sprintf("Samples with at least one gene: %d\n", nrow(gene_matrix_nonzero)))
  
  if (nrow(gene_matrix_nonzero) < 3) {
    cat("Not enough samples for PCoA analysis\n\n")
    return(NULL)
  }
  
  # Remove columns with all zeros
  col_sums <- colSums(gene_matrix_nonzero)
  gene_matrix_clean <- gene_matrix_nonzero[, col_sums > 0, drop = FALSE]
  
  cat(sprintf("Gene types present: %d\n", ncol(gene_matrix_clean)))
  
  # Convert to presence/absence for Jaccard, or use counts
  # Store original clean matrix for count mode
  if (args$use_counts) {
    cat("Using count numbers for distance calculation (Bray-Curtis distance)...\n")
    gene_matrix_for_dist <- gene_matrix_clean
    dist_matrix <- vegdist(gene_matrix_for_dist, method = "bray")
    # Also create binary version for saving (for compatibility with step2)
    gene_matrix_binary <- (gene_matrix_clean > 0) * 1
  } else {
    cat("Converting to presence/absence for Jaccard distance...\n")
    gene_matrix_binary <- (gene_matrix_clean > 0) * 1
    gene_matrix_for_dist <- gene_matrix_binary
    dist_matrix <- vegdist(gene_matrix_for_dist, method = "jaccard", binary = TRUE)
  }
  
  cat("Running PCoA...\n")
  set.seed(123)
  
  pcoa_result <- tryCatch({
    pcoa(dist_matrix, correction = "cailliez")
  }, error = function(e) {
    cat(sprintf("Error in PCoA: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(pcoa_result)) return(NULL)
  
  # Calculate variance explained
  eigenvalues <- pcoa_result$values$Eigenvalues
  eigenvalues[eigenvalues < 0] <- 0
  variance_explained <- eigenvalues / sum(eigenvalues) * 100
  
  cat(sprintf("PCoA Axis 1 variance explained: %.2f%%\n", variance_explained[1]))
  cat(sprintf("PCoA Axis 2 variance explained: %.2f%%\n", variance_explained[2]))
  
  # Extract scores
  pcoa_scores <- as.data.frame(pcoa_result$vectors[, 1:2])
  colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
  pcoa_scores$Sample_ID <- rownames(pcoa_scores)
  
  # Merge with sample info
  plot_data <- pcoa_scores %>%
    left_join(sample_info_nonzero, by = "Sample_ID")
  
  # Remove outliers using IQR method
  cat("Checking for outliers...\n")
  
  q1_x <- quantile(plot_data$PCoA1, 0.25, na.rm = TRUE)
  q3_x <- quantile(plot_data$PCoA1, 0.75, na.rm = TRUE)
  iqr_x <- q3_x - q1_x
  lower_x <- q1_x - 1.5 * iqr_x
  upper_x <- q3_x + 1.5 * iqr_x
  
  q1_y <- quantile(plot_data$PCoA2, 0.25, na.rm = TRUE)
  q3_y <- quantile(plot_data$PCoA2, 0.75, na.rm = TRUE)
  iqr_y <- q3_y - q1_y
  lower_y <- q1_y - 1.5 * iqr_y
  upper_y <- q3_y + 1.5 * iqr_y
  
  outlier_idx <- plot_data$PCoA1 < lower_x | plot_data$PCoA1 > upper_x |
                 plot_data$PCoA2 < lower_y | plot_data$PCoA2 > upper_y
  
  n_outliers <- sum(outlier_idx, na.rm = TRUE)
  max_outlier_pct <- if (args$type == "antidefense") 0.3 else 1.0
  
  if (n_outliers > 0 && n_outliers < nrow(plot_data) * max_outlier_pct) {
    cat(sprintf("Removing %d outliers (%.2f%%)\n", n_outliers, n_outliers/nrow(plot_data)*100))
    
    plot_data <- plot_data[!outlier_idx, ]
    keep_samples <- plot_data$Sample_ID
    
    # Update matrices based on whether we're using counts or binary
    if (args$use_counts) {
      # Filter the original clean matrix and recalculate binary
      gene_matrix_clean <- gene_matrix_clean[rownames(gene_matrix_clean) %in% keep_samples, , drop = FALSE]
      gene_matrix_clean <- gene_matrix_clean[, colSums(gene_matrix_clean) > 0, drop = FALSE]
      gene_matrix_binary <- (gene_matrix_clean > 0) * 1
      gene_matrix_for_dist <- gene_matrix_clean
    } else {
      gene_matrix_binary <- gene_matrix_binary[rownames(gene_matrix_binary) %in% keep_samples, , drop = FALSE]
      gene_matrix_binary <- gene_matrix_binary[, colSums(gene_matrix_binary) > 0, drop = FALSE]
      gene_matrix_for_dist <- gene_matrix_binary
    }
    
    # Check if we still have enough samples and gene types after outlier removal
    if (nrow(gene_matrix_for_dist) < 3) {
      cat("Not enough samples after outlier removal\n")
      return(NULL)
    }
    if (ncol(gene_matrix_for_dist) < 2) {
      cat("Not enough gene types after outlier removal\n")
      return(NULL)
    }
    
    cat(sprintf("Remaining samples: %d\n", nrow(plot_data)))
    
    # Re-run PCoA
    cat("Re-running PCoA on cleaned data...\n")
    if (args$use_counts) {
      dist_matrix <- vegdist(gene_matrix_for_dist, method = "bray")
    } else {
      dist_matrix <- vegdist(gene_matrix_for_dist, method = "jaccard", binary = TRUE)
    }
    
    pcoa_result_new <- tryCatch({
      pcoa(dist_matrix, correction = "cailliez")
    }, error = function(e) {
      if (args$type == "antidefense") {
        cat(sprintf("Warning: PCoA failed after outlier removal (%s)\n", e$message))
        cat("Keeping original PCoA results without outlier removal.\n")
        return(NULL)
      } else {
        stop(e)
      }
    })
    
    if (!is.null(pcoa_result_new)) {
      pcoa_result <- pcoa_result_new
      
      eigenvalues <- pcoa_result$values$Eigenvalues
      eigenvalues[eigenvalues < 0] <- 0
      variance_explained <- eigenvalues / sum(eigenvalues) * 100
      
      cat(sprintf("PCoA Axis 1 (after outlier removal): %.2f%%\n", variance_explained[1]))
      cat(sprintf("PCoA Axis 2 (after outlier removal): %.2f%%\n", variance_explained[2]))
      
      # Update scores
      pcoa_scores <- as.data.frame(pcoa_result$vectors[, 1:2])
      colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
      pcoa_scores$Sample_ID <- rownames(pcoa_scores)
      
      plot_data <- pcoa_scores %>%
        left_join(sample_info_nonzero, by = "Sample_ID")
    }
  } else if (n_outliers > 0) {
    cat(sprintf("Skipping outlier removal (too many outliers: %d, %.2f%%)\n", 
                n_outliers, n_outliers/nrow(plot_data)*100))
  } else {
    cat("No outliers detected.\n")
  }
  
  # Save PCoA coordinates with metadata
  coord_file <- file.path(output_dir, sprintf("%s_%s_coordinates.csv", output_prefix, contig_type))
  write.csv(plot_data, coord_file, row.names = FALSE)
  cat(sprintf("Coordinates saved to: %s\n", coord_file))
  
  # Save variance explained
  var_df <- data.frame(
    Axis = paste0("PCoA", 1:min(10, length(variance_explained))),
    Eigenvalue = eigenvalues[1:min(10, length(eigenvalues))],
    Variance_Explained_Pct = variance_explained[1:min(10, length(variance_explained))],
    Cumulative_Variance_Pct = cumsum(variance_explained)[1:min(10, length(variance_explained))]
  )
  var_file <- file.path(output_dir, sprintf("%s_%s_variance.csv", output_prefix, contig_type))
  write.csv(var_df, var_file, row.names = FALSE)
  cat(sprintf("Variance explained saved to: %s\n", var_file))
  
  # Save gene matrix (binary) for envfit in step2
  # Note: Always save binary version for compatibility with step2 plotting script
  matrix_file <- file.path(output_dir, sprintf("%s_%s_%s.csv", output_prefix, contig_type, matrix_suffix))
  gene_df <- as.data.frame(gene_matrix_binary)
  gene_df$Sample_ID <- rownames(gene_matrix_binary)
  gene_df <- gene_df %>% select(Sample_ID, everything())
  write.csv(gene_df, matrix_file, row.names = FALSE)
  if (args$use_counts) {
    cat(sprintf("Gene matrix (binary) saved to: %s (Note: PCoA was calculated using counts)\n\n", matrix_file))
  } else {
    cat(sprintf("Gene matrix (binary) saved to: %s\n\n", matrix_file))
  }
  
  return(list(coordinates = plot_data, variance = variance_explained))
}

# Function to process a single group or all data
process_group <- function(data_subset, group_name = NULL) {
  if (!is.null(group_name)) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Processing group: %s\n", group_name))
    cat("========================================\n\n")
    output_prefix <- sprintf("%s_%s", args$output, gsub("[^A-Za-z0-9_]", "_", group_name))
  } else {
    output_prefix <- args$output
  }
  
  result_chr <- run_pcoa_calculation(data_subset, "Chromosome", output_prefix, 
                                      existing_gene_cols, config$matrix_suffix, input_dir)
  result_pla <- run_pcoa_calculation(data_subset, "Plasmid", output_prefix, 
                                      existing_gene_cols, config$matrix_suffix, input_dir)
  
  return(list(chr = result_chr, pla = result_pla, prefix = output_prefix))
}

# Run PCoA calculation
cat("\n========================================\n")
cat(sprintf("Starting PCoA Calculation for %s\n", config$name))
if (!is.null(args$group)) {
  cat(sprintf("Grouping by: %s\n", args$group))
}
cat("========================================\n\n")

all_results <- list()

if (!is.null(args$group) && !is.null(groups)) {
  # Process each group separately
  for (group_val in groups) {
    df_group <- df_filtered[df_filtered[[args$group]] == group_val, ]
    cat(sprintf("\nGroup '%s': %d rows\n", group_val, nrow(df_group)))
    
    group_results <- process_group(df_group, group_val)
    all_results[[as.character(group_val)]] <- group_results
  }
} else {
  # Process all data together
  results <- process_group(df_filtered, NULL)
  all_results[["all"]] <- results
}

# Summary
cat("\n========================================\n")
cat("PCoA Calculation Complete!\n")
cat("========================================\n\n")

for (result_name in names(all_results)) {
  results <- all_results[[result_name]]
  if (!is.null(result_name) && result_name != "all") {
    cat(sprintf("Group: %s\n", result_name))
  }
  
  if (!is.null(results$chr)) {
    cat(sprintf("  Chromosome: %d samples, Axis1=%.2f%%, Axis2=%.2f%%\n", 
                nrow(results$chr$coordinates), results$chr$variance[1], results$chr$variance[2]))
  }
  if (!is.null(results$pla)) {
    cat(sprintf("  Plasmid: %d samples, Axis1=%.2f%%, Axis2=%.2f%%\n", 
                nrow(results$pla$coordinates), results$pla$variance[1], results$pla$variance[2]))
  }
  cat("\n")
}

cat("Output directory: ", input_dir, "\n")
cat("Output files saved to input file directory.\n")
cat(sprintf("\nRun 214_step2_pcoa_plot.R to create customized plots.\n"))

