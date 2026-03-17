#!/usr/bin/env Rscript

# =============================================================================
# Taxonomy Heatmap Generator
# Description:
#   Generate heatmaps of gene subtype distributions by selected taxonomy levels,
#   split by Contig_Type3 and faceted by Host.
# Usage:
#   Rscript 204_taxonomy_heatmap.R -i <input.csv> -o <output_dir> [-g <top_genes>] [--taxonomy-levels <levels>] [--top-kingdom <n>] [--top-phylum <n>] [--top-class <n>] [--top-order <n>] [--top-family <n>] [--top-genus-updated <n>] [--top-species <n>] [-w <width>] [-e <height>]
# Example:
#   Rscript 204_taxonomy_heatmap.R -i Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv -o Result/NCBI_4395_Batch/04_Gene_Taxonomy -g 50 --taxonomy-levels Kingdom_CRBC,Phylum_CRBC,Genus_CRBC_Updated --top-kingdom 20 --top-phylum 50 -w 15 -e 30
#   Rscript 204_taxonomy_heatmap.R -i Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv -o Result/NCBI_4395_Batch/04_Gene_Taxonomy/ -g 20 --taxonomy-levels Class_CRBC -w 12 -e 8
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(argparse)
})

# Step 0: Define taxonomy levels
taxonomy_levels_all <- c(
  "Kingdom_CRBC",
  "Phylum_CRBC",
  "Class_CRBC",
  "Order_CRBC",
  "Family_CRBC",
  "Genus_CRBC_Updated",
  "Species_CRBC"
)
contig_type3_order <- c("Chromosome", "Plasmid", "Virus")

# Step 1: Helper functions
parse_level_list <- function(levels_arg, available_levels) {
  levels <- strsplit(levels_arg, ",", fixed = TRUE)[[1]]
  levels <- trimws(levels)
  levels <- levels[levels != ""]
  levels <- unique(levels)

  invalid <- setdiff(levels, available_levels)
  if (length(invalid) > 0) {
    stop(sprintf(
      "Invalid taxonomy level(s): %s. Allowed levels: %s",
      paste(invalid, collapse = ", "),
      paste(available_levels, collapse = ", ")
    ))
  }

  if (length(levels) == 0) {
    stop("No valid taxonomy levels provided in --taxonomy-levels.")
  }

  levels
}

normalize_top_n <- function(value, arg_name) {
  if (is.null(value) || is.na(value) || value == -1) {
    return(NA_integer_)
  }
  if (value < 1) {
    stop(sprintf("%s must be >= 1 when provided.", arg_name))
  }
  as.integer(value)
}

# Step 1: Parse command line arguments
parser <- ArgumentParser(description = "Generate taxonomy heatmaps by Contig_Type3 and Host")
parser$add_argument(
  "-i", "--input", default = "Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv",
  help = "Input CSV file path (default: Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv)"
)
parser$add_argument(
  "-o", "--output", default = "Result/NCBI_4395_Batch/04_Gene_Taxonomy",
  help = "Output directory for PDF and CSV files (default: Result/NCBI_4395_Batch/04_Gene_Taxonomy)"
)
parser$add_argument(
  "-g", "--top-genes", type = "integer", default = 50,
  help = "Top N genes to keep per Contig_Type3 (default: 50)"
)
parser$add_argument(
  "--taxonomy-levels", default = paste(taxonomy_levels_all, collapse = ","),
  help = "Comma-separated taxonomy columns to use on X axis. Default is all 7 levels."
)
parser$add_argument(
  "--top-kingdom", type = "integer", default = -1,
  help = "Top N for Kingdom_CRBC (default: all)"
)
parser$add_argument(
  "--top-phylum", type = "integer", default = -1,
  help = "Top N for Phylum_CRBC (default: all)"
)
parser$add_argument(
  "--top-class", type = "integer", default = -1,
  help = "Top N for Class_CRBC (default: all)"
)
parser$add_argument(
  "--top-order", type = "integer", default = -1,
  help = "Top N for Order_CRBC (default: all)"
)
parser$add_argument(
  "--top-family", type = "integer", default = -1,
  help = "Top N for Family_CRBC (default: all)"
)
parser$add_argument(
  "--top-genus-updated", type = "integer", default = -1,
  help = "Top N for Genus_CRBC_Updated (default: all)"
)
parser$add_argument(
  "--top-species", type = "integer", default = -1,
  help = "Top N for Species_CRBC (default: all)"
)
parser$add_argument(
  "-w", "--fig-width", type = "double", default = 15,
  help = "Figure width in inches (default: 15)"
)
parser$add_argument(
  "-e", "--fig-height", type = "double", default = 30,
  help = "Figure height in inches (default: 30)"
)

args <- parser$parse_args()

if (args$top_genes < 1) {
  stop("--top-genes must be >= 1.")
}

if (args$fig_width <= 0 || args$fig_height <= 0) {
  stop("Both --fig-width and --fig-height must be > 0.")
}

args$top_kingdom <- normalize_top_n(args$top_kingdom, "--top-kingdom")
args$top_phylum <- normalize_top_n(args$top_phylum, "--top-phylum")
args$top_class <- normalize_top_n(args$top_class, "--top-class")
args$top_order <- normalize_top_n(args$top_order, "--top-order")
args$top_family <- normalize_top_n(args$top_family, "--top-family")
args$top_genus_updated <- normalize_top_n(args$top_genus_updated, "--top-genus-updated")
args$top_species <- normalize_top_n(args$top_species, "--top-species")

if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE, showWarnings = FALSE)
}

# Step 2: Resolve selected taxonomy levels and per-level TopN
selected_taxonomy_levels <- parse_level_list(args$taxonomy_levels, taxonomy_levels_all)
taxonomy_topn_map <- c(
  Kingdom_CRBC = args$top_kingdom,
  Phylum_CRBC = args$top_phylum,
  Class_CRBC = args$top_class,
  Order_CRBC = args$top_order,
  Family_CRBC = args$top_family,
  Genus_CRBC_Updated = args$top_genus_updated,
  Species_CRBC = args$top_species
)

# Step 3: Define defense filters
defense_exclude_list <- c("VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA")

# Step 4: Read and normalize input table
data <- read.csv(args$input, stringsAsFactors = FALSE, check.names = FALSE)
if (!("Contig_Type3" %in% names(data))) {
  stop("Missing required column: Contig_Type3")
}
data <- data %>%
  mutate(Contig_Type3 = if_else(Contig_Type3 %in% c("Virus", "Provirus"), "Virus", Contig_Type3))

unexpected_contig_type3 <- setdiff(unique(data$Contig_Type3), c(contig_type3_order, "", NA))
if (length(unexpected_contig_type3) > 0) {
  warning(sprintf(
    "Unexpected Contig_Type3 values found and removed: %s",
    paste(unexpected_contig_type3, collapse = ", ")
  ))
}

data <- data %>%
  filter(Contig_Type3 %in% contig_type3_order) %>%
  mutate(Contig_Type3 = factor(Contig_Type3, levels = contig_type3_order))

# Step 5: Validate required columns
required_gene_cols <- c("Defense_Subtype", "AntiDS_Type", "AMR_Type")
required_cols <- unique(c("Host", "Contig_Type3", required_gene_cols, selected_taxonomy_levels))
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "Missing required column(s): %s",
    paste(missing_cols, collapse = ", ")
  ))
}

# Step 6: Count genes by Host + Contig_Type3 + taxonomy
count_by_taxonomy <- function(df, gene_col, taxonomy_col, defense_filter = FALSE) {
  gene_df <- df %>%
    select(Host, Contig_Type3, taxonomy = all_of(taxonomy_col), all_of(gene_col)) %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    mutate(gene = strsplit(.data[[gene_col]], ",")) %>%
    unnest(gene) %>%
    mutate(gene = trimws(gene)) %>%
    filter(gene != "", !is.na(taxonomy), taxonomy != "")

  if (defense_filter) {
    gene_df <- gene_df %>%
      filter(!(gene %in% defense_exclude_list)) %>%
      filter(!str_detect(tolower(gene), "other"))
  }

  gene_df %>%
    count(Host, Contig_Type3, taxonomy, gene, name = "count")
}

# Step 7: Check whether a gene column should use defense filters
is_defense_like <- function(gene_col) {
  gene_col == "Defense_Subtype"
}

# Step 8: Draw heatmaps per Contig_Type3 and taxonomy level
create_heatmaps_by_contigtype <- function(df, gene_col, gene_name, taxonomy_col, top_genes, top_taxonomy, output_dir, fig_width, fig_height) {
  cat(sprintf("Processing %s with taxonomy %s...\n", gene_name, taxonomy_col))
  gene_counts <- count_by_taxonomy(
    df = df,
    gene_col = gene_col,
    taxonomy_col = taxonomy_col,
    defense_filter = is_defense_like(gene_col)
  )
  contig_types <- intersect(contig_type3_order, as.character(unique(gene_counts$Contig_Type3)))

  for (contig_type in contig_types) {
    cat(sprintf("  Plotting for Contig_Type3: %s, taxonomy: %s\n", contig_type, taxonomy_col))
    subdata <- gene_counts %>% filter(Contig_Type3 == contig_type)

    if (nrow(subdata) == 0) {
      cat(sprintf("    No data for Contig_Type3 %s, skipping...\n", contig_type))
      next
    }

    top_gene_values <- subdata %>%
      group_by(gene) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = top_genes) %>%
      pull(gene)

    plot_data <- subdata %>%
      filter(gene %in% top_gene_values)

    if (!is.na(top_taxonomy)) {
      top_taxonomy_values <- plot_data %>%
        group_by(taxonomy) %>%
        summarise(total = sum(count), .groups = "drop") %>%
        arrange(desc(total)) %>%
        slice_head(n = top_taxonomy) %>%
        pull(taxonomy)

      plot_data <- plot_data %>%
        filter(taxonomy %in% top_taxonomy_values)
    }

    if (nrow(plot_data) == 0) {
      cat(sprintf("    No rows left after top filtering for %s, skipping...\n", contig_type))
      next
    }

    gene_levels <- subdata %>%
      group_by(gene) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      filter(gene %in% top_gene_values) %>%
      arrange(desc(total)) %>%
      pull(gene)
    plot_data$gene <- factor(plot_data$gene, levels = gene_levels)

    taxonomy_levels <- plot_data %>%
      group_by(taxonomy) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      arrange(desc(total)) %>%
      pull(taxonomy)
    plot_data$taxonomy <- factor(plot_data$taxonomy, levels = taxonomy_levels)

    n_host <- length(unique(plot_data$Host))
    facet_ncol <- 4
    facet_nrow <- max(1, ceiling(n_host / facet_ncol))
    safe_contig <- gsub("[^A-Za-z0-9]", "_", contig_type)
    safe_taxonomy <- gsub("[^A-Za-z0-9]", "_", taxonomy_col)
    taxonomy_label <- ifelse(is.na(top_taxonomy), "ALL", paste0("TOP", top_taxonomy))

    pdf_name <- sprintf(
      "%s_%s_%s_heatmap_by_host_TOP%d_%s.pdf",
      gene_name, safe_taxonomy, safe_contig, top_genes, taxonomy_label
    )
    pdf_path <- file.path(output_dir, pdf_name)

    p <- ggplot(plot_data, aes(x = taxonomy, y = gene, fill = count)) +
      geom_tile(color = "white", linewidth = 0.12) +
      scale_fill_gradient(
        low = "white", high = "#6566aa",
        name = "Count",
        trans = "log1p"
      ) +
      facet_wrap(~Host, ncol = facet_ncol, nrow = facet_nrow, scales = "fixed") +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.position = "right",
        panel.grid = element_blank()
      ) +
      labs(
        title = sprintf(
          "%s Top%d by %s (%s, Contig_Type3: %s)",
          gene_name, top_genes, taxonomy_col, taxonomy_label, contig_type
        ),
        x = taxonomy_col,
        y = gene_name
      )

    ggsave(pdf_path, plot = p, width = fig_width, height = fig_height, limitsize = FALSE)

    csv_name <- sprintf(
      "%s_%s_%s_counts_by_host_TOP%d_%s.csv",
      gene_name, safe_taxonomy, safe_contig, top_genes, taxonomy_label
    )
    csv_path <- file.path(output_dir, csv_name)
    write_csv(plot_data, csv_path)
  }
}

# Step 9: Run all selected taxonomy levels and gene systems
gene_systems <- list(
  list(col = "Defense_Subtype", name = "Defense_Subtype"),
  list(col = "AntiDS_Type", name = "AntiDS_Type"),
  list(col = "AMR_Type", name = "AMR_Type")
)

for (taxonomy_col in selected_taxonomy_levels) {
  current_top_taxonomy <- taxonomy_topn_map[[taxonomy_col]]
  for (gene_sys in gene_systems) {
    create_heatmaps_by_contigtype(
      df = data,
      gene_col = gene_sys$col,
      gene_name = gene_sys$name,
      taxonomy_col = taxonomy_col,
      top_genes = args$top_genes,
      top_taxonomy = current_top_taxonomy,
      output_dir = args$output,
      fig_width = args$fig_width,
      fig_height = args$fig_height
    )
  }
}

cat("\nAll heatmaps completed!\n")
