#!/usr/bin/env Rscript
# Remove parallel edges (duplicate edges between the same pair of nodes)
# For undirected graphs, edges (A,B) and (B,A) are considered the same
# Keep the edge with the highest jaccard value if multiple edges exist

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript remove_parallel_edges.R input.tsv [output.tsv]")
}

input_file <- args[1]
if (length(args) >= 2) {
  output_file <- args[2]
} else {
  # Default output: add _no_parallel suffix
  output_file <- sub("\\.tsv$", "_no_parallel.tsv", input_file)
  if (output_file == input_file) {
    output_file <- paste0(input_file, "_no_parallel.tsv")
  }
}

cat(sprintf("Reading file: %s\n", input_file))
edges_df <- read.delim(input_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Original edges: %d\n", nrow(edges_df)))

# Check required columns
if (!("g1" %in% colnames(edges_df) && "g2" %in% colnames(edges_df))) {
  stop("Edge file must have 'g1' and 'g2' columns")
}

# Normalize edge pairs: ensure g1 < g2 (for undirected graph)
# This makes (A,B) and (B,A) the same
normalize_edge_pair <- function(g1, g2) {
  # Create a data frame and sort each row
  pairs <- data.frame(g1 = g1, g2 = g2, stringsAsFactors = FALSE)
  for (i in 1:nrow(pairs)) {
    if (pairs$g1[i] > pairs$g2[i]) {
      # Swap
      temp <- pairs$g1[i]
      pairs$g1[i] <- pairs$g2[i]
      pairs$g2[i] <- temp
    }
  }
  return(pairs)
}

cat("Normalizing edge pairs (ensuring g1 < g2 for undirected graph)...\n")
normalized_pairs <- normalize_edge_pair(edges_df$g1, edges_df$g2)
edges_df$g1_norm <- normalized_pairs$g1
edges_df$g2_norm <- normalized_pairs$g2

# Create a unique edge identifier
edges_df$edge_id <- paste(edges_df$g1_norm, edges_df$g2_norm, sep = "|")

# Count parallel edges
edge_counts <- table(edges_df$edge_id)
parallel_edges <- sum(edge_counts > 1)
cat(sprintf("Found %d unique edge pairs\n", length(edge_counts)))
cat(sprintf("Found %d edge pairs with parallel edges\n", parallel_edges))

# Remove parallel edges: keep the one with highest jaccard value
if ("jaccard" %in% colnames(edges_df)) {
  cat("Removing parallel edges, keeping edge with highest jaccard value...\n")
  
  # Sort by edge_id and jaccard (descending)
  edges_df <- edges_df[order(edges_df$edge_id, -edges_df$jaccard), ]
  
  # Keep first occurrence of each edge_id (which has highest jaccard)
  edges_df_unique <- edges_df[!duplicated(edges_df$edge_id), ]
  
} else {
  cat("No jaccard column found, keeping first occurrence of each edge pair...\n")
  edges_df <- edges_df[order(edges_df$edge_id), ]
  edges_df_unique <- edges_df[!duplicated(edges_df$edge_id), ]
}

# Remove temporary columns
edges_df_unique$g1_norm <- NULL
edges_df_unique$g2_norm <- NULL
edges_df_unique$edge_id <- NULL

cat(sprintf("Edges after removing parallel edges: %d\n", nrow(edges_df_unique)))
cat(sprintf("Removed %d parallel edges\n", nrow(edges_df) - nrow(edges_df_unique)))

# Write output
cat(sprintf("Writing output to: %s\n", output_file))
write.table(edges_df_unique, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Parallel edges removal completed successfully!\n")
cat(sprintf("Output saved to: %s\n", output_file))
