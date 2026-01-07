#!/usr/bin/env Rscript

# =============================================================================
# Step 2: PCoA Plotting (Generic)
# Description: Create customized PCoA plots from pre-calculated coordinates
#              Supports: Defense, AMR, AntiDefense
#              Supports: color by taxonomy, shape by Host, gene arrows,
#              confidence ellipses, and more
# Usage: Rscript 214_step2_pcoa_plot.R -p <prefix> -t <type> [options]
#
# Arguments:
#   -p: Input file prefix (e.g., PCoA_Data, PCoA_AMR_Data, PCoA_AntiDS_Data)
#   -t: Gene type: 'defense', 'amr', or 'antidefense'
#   -c: Optional color by variable (default: Class_CRBC if not specified)
#   -s: Optional shape by variable (default: Host if not specified, use 'NULL' for none)
#   -e: Show ellipse (default: TRUE)
#   -a: Show arrows (default: TRUE)
#   -o: Output suffix (default: final)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Create customized PCoA plots")
parser$add_argument("-p", "--prefix", required = TRUE, 
                    help = "Input file prefix (e.g., PCoA_Data, PCoA_AMR_Data)")
parser$add_argument("-t", "--type", required = TRUE,
                    choices = c("defense", "amr", "antidefense"),
                    help = "Gene type: 'defense', 'amr', or 'antidefense'")
parser$add_argument("-c", "--color-by", default = NULL,
                    help = "Optional color by variable (default: Class_CRBC if not specified)")
parser$add_argument("-s", "--shape-by", default = NULL,
                    help = "Optional shape by variable (default: Host if not specified, use 'NULL' for none)")
parser$add_argument("-e", "--ellipse", action = "store_true", default = TRUE,
                    help = "Show confidence ellipses")
parser$add_argument("-a", "--arrows", action = "store_true", default = TRUE,
                    help = "Show gene arrows (envfit)")
parser$add_argument("-o", "--output-suffix", default = "final",
                    help = "Output file suffix (default: final)")

args <- parser$parse_args()

# Configuration
type_config <- list(
  defense = list(
    name = "Defense Systems",
    matrix_suffix = "defense_matrix",
    label_prefix = "Defense_"
  ),
  amr = list(
    name = "AMR Types",
    matrix_suffix = "amr_matrix",
    label_prefix = "AMR_"
  ),
  antidefense = list(
    name = "Anti-Defense Systems",
    matrix_suffix = "antids_matrix",
    label_prefix = "AntiDS_"
  )
)

config <- type_config[[args$type]]

# Settings - use defaults if not provided
COLOR_BY <- if (is.null(args$color_by)) "Class_CRBC" else args$color_by
SHAPE_BY <- if (is.null(args$shape_by)) "Host" else if (args$shape_by == "NULL") NULL else args$shape_by
SHOW_ELLIPSE <- args$ellipse
SHOW_ARROWS <- args$arrows
ELLIPSE_LEVEL <- 0.95
N_TOP_ARROWS <- 15
ARROW_PVAL <- 0.05
ARROW_R2 <- 0.05
POINT_SIZE <- 2.5
POINT_ALPHA <- 0.6
OUTPUT_SUFFIX <- args$output_suffix

# Color palettes
class_colors <- c(
  "Actinomycetia" = "#98df8a",
  "Alphaproteobacteria" = "#aec7e8",
  "Bacilli" = "#ff7f0e",
  "Thermoleophilia" = "#ff9896",
  "Bacteroidia" = "#d62728",
  "Gammaproteobacteria" = "#ffbb78",
  "Deinococci" = "#1f77b4",
  "Acidimicrobiia" = "#2ca02c",
  "Campylobacteria" = "#9467bd"
)

host_shapes <- c(
  "Alfalfa" = 16,
  "Rice" = 17,
  "Wheat" = 15,
  "Maize" = 18
)

# Function to calculate envfit arrows
calculate_envfit_arrows <- function(plot_data, gene_matrix, label_prefix) {
  cat(sprintf("Calculating %s arrows (envfit)...\n", config$name))
  
  common_samples <- intersect(plot_data$Sample_ID, gene_matrix$Sample_ID)
  plot_data_matched <- plot_data %>% filter(Sample_ID %in% common_samples)
  gene_matched <- gene_matrix %>% filter(Sample_ID %in% common_samples)
  
  plot_data_matched <- plot_data_matched[match(common_samples, plot_data_matched$Sample_ID), ]
  gene_matched <- gene_matched[match(common_samples, gene_matched$Sample_ID), ]
  
  pcoa_coords <- as.matrix(plot_data_matched[, c("PCoA1", "PCoA2")])
  gene_cols <- setdiff(colnames(gene_matched), "Sample_ID")
  gene_mat <- as.matrix(gene_matched[, gene_cols])
  
  min_presence <- ceiling(nrow(gene_mat) * 0.05)
  col_sums <- colSums(gene_mat)
  common_genes <- names(col_sums[col_sums >= min_presence])
  
  cat(sprintf("%s types present in >=5%% samples: %d\n", config$name, length(common_genes)))
  
  if (length(common_genes) < 2) {
    cat(sprintf("Not enough common %s types for arrows\n", config$name))
    return(data.frame())
  }
  
  gene_mat <- gene_mat[, common_genes, drop = FALSE]
  envfit_df <- data.frame()
  
  for (gene in colnames(gene_mat)) {
    gene_vec <- gene_mat[, gene]
    if (sd(gene_vec) == 0) next
    
    cor1 <- cor(gene_vec, pcoa_coords[, 1])
    cor2 <- cor(gene_vec, pcoa_coords[, 2])
    r2 <- cor1^2 + cor2^2
    
    n_perm <- 999
    perm_r2 <- numeric(n_perm)
    for (i in 1:n_perm) {
      perm_vec <- sample(gene_vec)
      perm_cor1 <- cor(perm_vec, pcoa_coords[, 1])
      perm_cor2 <- cor(perm_vec, pcoa_coords[, 2])
      perm_r2[i] <- perm_cor1^2 + perm_cor2^2
    }
    pval <- (sum(perm_r2 >= r2) + 1) / (n_perm + 1)
    
    vec_length <- sqrt(cor1^2 + cor2^2)
    if (vec_length > 0) {
      norm_cor1 <- cor1 / vec_length * sqrt(r2)
      norm_cor2 <- cor2 / vec_length * sqrt(r2)
    } else {
      norm_cor1 <- 0
      norm_cor2 <- 0
    }
    
    envfit_df <- rbind(envfit_df, data.frame(
      gene = gene,
      PCoA1 = norm_cor1,
      PCoA2 = norm_cor2,
      pval = pval,
      r2 = r2
    ))
  }
  
  sig_vectors <- envfit_df %>%
    filter(pval < ARROW_PVAL & r2 > ARROW_R2) %>%
    arrange(desc(r2))
  
  cat(sprintf("Significant vectors (p<%.2f, r2>%.2f): %d\n", 
              ARROW_PVAL, ARROW_R2, nrow(sig_vectors)))
  
  return(sig_vectors)
}

# Function to generate color palette
generate_colors <- function(categories) {
  n <- length(categories)
  
  if (n <= 9 && all(categories %in% names(class_colors))) {
    return(class_colors[categories])
  }
  
  if (n <= 12) {
    colors <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", 
                "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
                "#ccebc5", "#ffed6f")
  } else {
    colors <- c(
      "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
      "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f",
      "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
      "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928",
      "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628",
      "#f781bf", "#999999", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
      "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"
    )
  }
  
  return(setNames(colors[1:n], categories))
}

# Main plotting function
create_pcoa_plot <- function(contig_type, input_dir, file_prefix, base_prefix) {
  cat(sprintf("\n=== Creating plot for %s ===\n", contig_type))
  
  coord_file <- file.path(input_dir, sprintf("%s_%s_coordinates.csv", file_prefix, contig_type))
  var_file <- file.path(input_dir, sprintf("%s_%s_variance.csv", file_prefix, contig_type))
  matrix_file <- file.path(input_dir, sprintf("%s_%s_%s.csv", file_prefix, contig_type, config$matrix_suffix))
  
  if (!file.exists(coord_file)) {
    cat(sprintf("Coordinate file not found: %s\n", coord_file))
    cat("Please run 213_step1_pcoa_calculate.R first.\n")
    return(NULL)
  }
  
  plot_data <- read.csv(coord_file, stringsAsFactors = FALSE)
  cat(sprintf("Loaded %d samples\n", nrow(plot_data)))
  
  var_axis1 <- NA
  var_axis2 <- NA
  if (file.exists(var_file)) {
    var_data <- read.csv(var_file, stringsAsFactors = FALSE)
    var_axis1 <- var_data$Variance_Explained_Pct[1]
    var_axis2 <- var_data$Variance_Explained_Pct[2]
  }
  cat(sprintf("Variance: Axis1=%.2f%%, Axis2=%.2f%%\n", var_axis1, var_axis2))
  
  top_vectors <- data.frame()
  if (SHOW_ARROWS && file.exists(matrix_file)) {
    gene_matrix <- read.csv(matrix_file, stringsAsFactors = FALSE)
    sig_vectors <- calculate_envfit_arrows(plot_data, gene_matrix, config$label_prefix)
    top_vectors <- head(sig_vectors, N_TOP_ARROWS)
    
    if (nrow(sig_vectors) > 0) {
      # Extract group name from file_prefix if it exists
      if (file_prefix != base_prefix && grepl("_", file_prefix, fixed = TRUE)) {
        group_name <- sub(paste0("^", base_prefix, "_"), "", file_prefix)
        envfit_file <- file.path(input_dir, sprintf("PCoA_%s_%s_%s_%s_envfit.csv", 
                                toupper(substr(args$type, 1, 1)), 
                                group_name, contig_type, OUTPUT_SUFFIX))
      } else {
        envfit_file <- file.path(input_dir, sprintf("PCoA_%s_%s_%s_envfit.csv", 
                                toupper(substr(args$type, 1, 1)), 
                                contig_type, OUTPUT_SUFFIX))
      }
      write.csv(sig_vectors, envfit_file, row.names = FALSE)
      cat(sprintf("Envfit results saved to: %s\n", envfit_file))
    }
  }
  
  color_var <- COLOR_BY
  all_categories <- sort(unique(plot_data[[color_var]][!is.na(plot_data[[color_var]])]))
  color_palette <- generate_colors(all_categories)
  
  x_range <- range(plot_data$PCoA1)
  y_range <- range(plot_data$PCoA2)
  x_pad <- diff(x_range) * 0.15
  y_pad <- diff(y_range) * 0.15
  x_limits <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  y_limits <- c(y_range[1] - y_pad, y_range[2] + y_pad)
  
  arrow_scale <- min(diff(x_range), diff(y_range)) * 0.5
  
  if (nrow(top_vectors) > 0) {
    top_vectors$x_end <- top_vectors$PCoA1 * arrow_scale
    top_vectors$y_end <- top_vectors$PCoA2 * arrow_scale
    top_vectors$gene_label <- gsub(config$label_prefix, "", top_vectors$gene)
  }
  
  cat("Building plot...\n")
  
  if (!is.null(SHAPE_BY) && SHAPE_BY %in% colnames(plot_data)) {
    p <- ggplot(plot_data, aes_string(x = "PCoA1", y = "PCoA2", 
                                       color = color_var, shape = SHAPE_BY))
  } else {
    p <- ggplot(plot_data, aes_string(x = "PCoA1", y = "PCoA2", color = color_var))
  }
  
  if (SHOW_ELLIPSE) {
    p <- p + stat_ellipse(aes_string(fill = color_var, group = color_var), 
                          geom = "polygon", 
                          alpha = 0.12, 
                          level = ELLIPSE_LEVEL,
                          type = "t",
                          show.legend = FALSE)
  }
  
  p <- p + geom_point(alpha = POINT_ALPHA, size = POINT_SIZE)
  p <- p + scale_color_manual(values = color_palette, name = gsub("_CRBC", "", color_var))
  p <- p + scale_fill_manual(values = color_palette)
  
  if (!is.null(SHAPE_BY) && SHAPE_BY == "Host") {
    p <- p + scale_shape_manual(values = host_shapes, name = "Host")
  }
  
  p <- p + labs(
    title = sprintf("PCoA of %s - %s", config$name, contig_type),
    subtitle = sprintf("Jaccard distance | Axis1: %.2f%%, Axis2: %.2f%% | n = %d samples", 
                       var_axis1, var_axis2, nrow(plot_data)),
    x = sprintf("PCoA1 (%.2f%%)", var_axis1),
    y = sprintf("PCoA2 (%.2f%%)", var_axis2)
  )
  
  p <- p + theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits)
  
  p <- p + guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2, override.aes = list(size = 3)),
    fill = "none"
  )
  
  if (nrow(top_vectors) > 0) {
    p <- p +
      geom_segment(data = top_vectors,
                   aes(x = 0, y = 0, xend = x_end, yend = y_end),
                   arrow = arrow(length = unit(0.2, "cm")),
                   color = "darkgray", linewidth = 0.8, inherit.aes = FALSE) +
      geom_text_repel(data = top_vectors,
                      aes(x = x_end, y = y_end, label = gene_label),
                      color = "black", size = 3, inherit.aes = FALSE,
                      max.overlaps = 20, segment.color = "gray70",
                      fontface = "italic")
  }
  
  # Extract group name from file_prefix if it exists
  # file_prefix format: {base_prefix} or {base_prefix}_{group_name}
  if (file_prefix != base_prefix && grepl("_", file_prefix, fixed = TRUE)) {
    # Extract group name by removing base_prefix and the underscore
    group_name <- sub(paste0("^", base_prefix, "_"), "", file_prefix)
    output_file <- file.path(input_dir, sprintf("PCoA_%s_%s_%s_%s.pdf", 
                           toupper(substr(args$type, 1, 1)), 
                           group_name, contig_type, OUTPUT_SUFFIX))
  } else {
    output_file <- file.path(input_dir, sprintf("PCoA_%s_%s_%s.pdf", 
                           toupper(substr(args$type, 1, 1)), 
                           contig_type, OUTPUT_SUFFIX))
  }
  
  ggsave(output_file, p, width = 12, height = 10)
  cat(sprintf("Plot saved to: %s\n", output_file))
  
  return(p)
}

# Function to find all coordinate files matching the prefix pattern
find_coordinate_files <- function(prefix, input_dir) {
  # Pattern: {prefix} or {prefix}_{group}_Chromosome/Plasmid_coordinates.csv
  pattern <- sprintf("^%s(_[^_]+)?_(Chromosome|Plasmid)_coordinates\\.csv$", prefix)
  
  all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = FALSE)
  coord_files <- all_files[grepl(pattern, all_files)]
  
  if (length(coord_files) == 0) {
    return(character(0))
  }
  
  # Extract unique prefixes (base + optional group)
  file_prefixes <- unique(sub("_(Chromosome|Plasmid)_coordinates\\.csv$", "", coord_files))
  
  return(file_prefixes)
}

# Determine input directory from prefix
# Check if prefix contains a path (directory separator)
if (grepl("[/\\\\]", args$prefix)) {
  # Prefix contains path, extract directory and prefix
  input_dir <- dirname(args$prefix)
  base_prefix <- basename(args$prefix)
} else {
  # Prefix is just a name, use current directory
  input_dir <- "."
  base_prefix <- args$prefix
}

# Normalize input directory
if (!dir.exists(input_dir)) {
  stop(sprintf("Directory does not exist: %s\n", input_dir))
}
input_dir <- normalizePath(input_dir, mustWork = TRUE)
cat(sprintf("Input directory: %s\n", input_dir))
cat(sprintf("Base prefix: %s\n", base_prefix))

# Find all coordinate files matching the prefix pattern
file_prefixes <- find_coordinate_files(base_prefix, input_dir)

if (length(file_prefixes) == 0) {
  stop(sprintf("Cannot find any coordinate files with prefix '%s' in directory '%s'.\n", base_prefix, input_dir))
}

cat(sprintf("Found %d file prefix(es): %s\n", length(file_prefixes), paste(file_prefixes, collapse = ", ")))

# Run plotting
cat("\n========================================\n")
cat("Step 2: PCoA Plotting\n")
cat("========================================\n")
cat(sprintf("\nSettings:\n"))
cat(sprintf("  Type: %s\n", config$name))
cat(sprintf("  Color by: %s\n", COLOR_BY))
cat(sprintf("  Shape by: %s\n", ifelse(is.null(SHAPE_BY), "None", SHAPE_BY)))
cat(sprintf("  Show ellipse: %s (level=%.0f%%)\n", SHOW_ELLIPSE, ELLIPSE_LEVEL*100))
cat(sprintf("  Show arrows: %s\n", SHOW_ARROWS))

# Process each file prefix
all_plots <- list()
for (file_prefix in file_prefixes) {
  cat(sprintf("\n>>> Processing prefix: %s <<<\n", file_prefix))
  
  p_chr <- create_pcoa_plot("Chromosome", input_dir, file_prefix, base_prefix)
  p_pla <- create_pcoa_plot("Plasmid", input_dir, file_prefix, base_prefix)
  
  all_plots[[file_prefix]] <- list(chr = p_chr, pla = p_pla)
}

cat("\n========================================\n")
cat("Plotting Complete!\n")
cat("========================================\n")
cat(sprintf("\nOutput directory: %s\n", input_dir))
cat(sprintf("Output files:\n"))
for (file_prefix in file_prefixes) {
  if (file_prefix != base_prefix) {
    group_name <- sub(paste0("^", base_prefix, "_"), "", file_prefix)
    cat(sprintf("  Group '%s':\n", group_name))
    cat(sprintf("    - PCoA_%s_%s_Chromosome_%s.pdf\n", toupper(substr(args$type, 1, 1)), group_name, OUTPUT_SUFFIX))
    cat(sprintf("    - PCoA_%s_%s_Plasmid_%s.pdf\n", toupper(substr(args$type, 1, 1)), group_name, OUTPUT_SUFFIX))
  } else {
    cat(sprintf("  - PCoA_%s_Chromosome_%s.pdf\n", toupper(substr(args$type, 1, 1)), OUTPUT_SUFFIX))
    cat(sprintf("  - PCoA_%s_Plasmid_%s.pdf\n", toupper(substr(args$type, 1, 1)), OUTPUT_SUFFIX))
  }
}
cat("\nTo customize the plot, use command line arguments.\n")

