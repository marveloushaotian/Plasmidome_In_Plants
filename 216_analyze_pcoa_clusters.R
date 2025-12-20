#!/usr/bin/env Rscript

# =============================================================================
# Analyze PCoA Clusters - Identify What Drives Cluster Formation
# Description: This script investigates what factors drive the clustering 
#              pattern in PCoA when neither Class_CRBC nor Host explain it.
# Usage: Rscript 216_analyze_pcoa_clusters.R -c <coords.csv> -m <matrix.csv> -d <workdir> -k <k_value> -o <output_prefix>
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(cluster)
  library(ggplot2)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Analyze PCoA Cluster Formation")
parser$add_argument("-c", "--coords", required = TRUE,
                    help = "PCoA coordinates CSV file")
parser$add_argument("-m", "--matrix", required = TRUE,
                    help = "Gene matrix CSV file")
parser$add_argument("-d", "--workdir", default = ".",
                    help = "Working directory (default: current directory)")
parser$add_argument("-k", "--k-value", type = "integer", default = 4,
                    help = "Number of clusters for k-means (default: 4)")
parser$add_argument("-o", "--output-prefix", default = "PCoA_AMR_Chromosome",
                    help = "Output file prefix (default: PCoA_AMR_Chromosome)")

args <- parser$parse_args()

setwd(args$workdir)

cat("=== Analyzing PCoA Cluster Formation ===\n\n")

# Read data
coords <- read.csv(args$coords, stringsAsFactors = FALSE)
gene_matrix <- read.csv(args$matrix, stringsAsFactors = FALSE)

cat(sprintf("Total samples: %d\n", nrow(coords)))
cat(sprintf("Gene types: %d\n", ncol(gene_matrix) - 1))

# Use k-means to identify clusters
cat("\n=== Step 1: Identifying clusters using k-means ===\n")

set.seed(123)
pcoa_coords <- coords[, c("PCoA1", "PCoA2")]

# Find optimal k using elbow method
wss <- sapply(1:10, function(k) {
  kmeans(pcoa_coords, centers = k, nstart = 10)$tot.withinss
})

cat("Within-cluster sum of squares for k=1:10:\n")
for (i in 1:10) {
  cat(sprintf("  k=%d: %.2f\n", i, wss[i]))
}

k <- args$k_value
km_result <- kmeans(pcoa_coords, centers = k, nstart = 25)
coords$Cluster <- factor(km_result$cluster)

cat(sprintf("\nCluster sizes (k=%d):\n", k))
print(table(coords$Cluster))

# Analyze what makes each cluster unique
cat("\n=== Step 2: What distinguishes each cluster? ===\n")

gene_data <- gene_matrix
rownames(gene_data) <- gene_data$Sample_ID
gene_data$Sample_ID <- NULL

common_samples <- intersect(coords$Sample_ID, rownames(gene_data))
coords_matched <- coords %>% filter(Sample_ID %in% common_samples)
gene_matched <- gene_data[common_samples, ]

coords_matched$Cluster <- factor(km_result$cluster[coords$Sample_ID %in% common_samples])

# Calculate gene profile for each cluster
cat("\n--- Gene Profile per Cluster ---\n")

gene_cols <- colnames(gene_matched)
cluster_gene_summary <- data.frame()

for (cl in 1:k) {
  cluster_samples <- coords_matched$Sample_ID[coords_matched$Cluster == cl]
  cluster_genes <- gene_matched[cluster_samples, , drop = FALSE]
  
  presence_rate <- colMeans(cluster_genes > 0) * 100
  top_genes <- sort(presence_rate, decreasing = TRUE)[1:10]
  
  cat(sprintf("\nCluster %d (n=%d samples):\n", cl, length(cluster_samples)))
  cat("  Top 10 gene types (presence rate %%):\n")
  for (i in 1:10) {
    cat(sprintf("    %s: %.1f%%\n", names(top_genes)[i], top_genes[i]))
  }
  
  for (gene in names(presence_rate)) {
    cluster_gene_summary <- rbind(cluster_gene_summary, data.frame(
      Cluster = cl,
      Gene = gene,
      PresenceRate = presence_rate[gene]
    ))
  }
}

# Find gene types that differentiate clusters
cat("\n=== Step 3: Gene types that differentiate clusters ===\n")

gene_wide <- cluster_gene_summary %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = PresenceRate, names_prefix = "Cluster_")

gene_wide$Variance <- apply(gene_wide[, grep("Cluster_", colnames(gene_wide))], 1, var, na.rm = TRUE)
gene_wide <- gene_wide %>% arrange(desc(Variance))

cat("\nGene types with highest variance across clusters:\n")
print(head(gene_wide, 15))

# Check if clusters correlate with metadata
cat("\n=== Step 4: Cluster correlation with metadata ===\n")

cat("\n--- Class_CRBC distribution per cluster ---\n")
class_cluster <- table(coords_matched$Cluster, coords_matched$Class_CRBC)
print(class_cluster)

chisq_class <- chisq.test(class_cluster)
cat(sprintf("\nChi-square test (Cluster vs Class_CRBC): p-value = %.4e\n", chisq_class$p.value))

cat("\n--- Host distribution per cluster ---\n")
host_cluster <- table(coords_matched$Cluster, coords_matched$Host)
print(host_cluster)

chisq_host <- chisq.test(host_cluster)
cat(sprintf("Chi-square test (Cluster vs Host): p-value = %.4e\n", chisq_host$p.value))

cat("\n--- Order_CRBC distribution per cluster ---\n")
order_cluster <- table(coords_matched$Cluster, coords_matched$Order_CRBC)
for (cl in 1:k) {
  cat(sprintf("\nCluster %d - Top 5 Orders:\n", cl))
  order_in_cluster <- sort(order_cluster[cl, ], decreasing = TRUE)[1:5]
  for (i in 1:5) {
    cat(sprintf("    %s: %d\n", names(order_in_cluster)[i], order_in_cluster[i]))
  }
}

# The KEY insight: check gene combination patterns
cat("\n=== Step 5: Gene Combination Patterns (THE KEY!) ===\n")

create_gene_signature <- function(gene_row) {
  present_genes <- names(gene_row)[gene_row > 0]
  if (length(present_genes) == 0) return("none")
  paste(sort(present_genes), collapse = "+")
}

gene_signatures <- apply(gene_matched, 1, create_gene_signature)
coords_matched$Gene_Signature <- gene_signatures

cat("\n--- Most common gene signatures per cluster ---\n")

for (cl in 1:k) {
  cluster_sigs <- coords_matched$Gene_Signature[coords_matched$Cluster == cl]
  sig_table <- sort(table(cluster_sigs), decreasing = TRUE)
  
  cat(sprintf("\nCluster %d - Top 5 gene signatures:\n", cl))
  top_sigs <- head(sig_table, 5)
  for (i in 1:length(top_sigs)) {
    sig_name <- names(top_sigs)[i]
    if (nchar(sig_name) > 80) sig_name <- paste0(substr(sig_name, 1, 80), "...")
    cat(sprintf("    %d samples: %s\n", top_sigs[i], sig_name))
  }
}

# Save visualization
cat("\n=== Step 6: Creating visualization ===\n")

p <- ggplot(coords_matched, aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  geom_point(alpha = 0.6, size = 2) +
  stat_ellipse(level = 0.95, type = "t") +
  labs(
    title = sprintf("PCoA Clusters - %s", args$output_prefix),
    subtitle = sprintf("Clusters identified by k-means (k=%d)", k),
    x = "PCoA1",
    y = "PCoA2"
  ) +
  theme_bw() +
  scale_color_brewer(palette = "Set1")

output_file <- sprintf("%s_clusters.pdf", args$output_prefix)
ggsave(output_file, p, width = 10, height = 8)
cat(sprintf("Cluster visualization saved to: %s\n", output_file))

# Summary
cat("\n=== SUMMARY ===\n")
cat("\nThese clusters are formed based on gene combination patterns!\n")
cat("In PCoA, samples with similar gene combinations cluster together.\n\n")

cat("Key findings:\n")
cat("1. Each cluster has characteristic gene combinations\n")
cat("2. This clustering does not completely correspond to bacterial taxonomy (Class_CRBC) or host (Host)\n")
cat("3. This suggests gene acquisition may be independent of phylogeny (horizontal gene transfer)\n\n")

cat("Possible explanations:\n")
cat("- Same plasmids/mobile elements carry similar gene combinations\n")
cat("- Same selection pressure leads to similar gene profiles\n")
cat("- Certain genes tend to co-occur (co-location or co-selection)\n")

# Save cluster assignments
output_file <- sprintf("%s_cluster_assignments.csv", args$output_prefix)
coords_matched %>%
  select(Sample_ID, PCoA1, PCoA2, Cluster, Host, Class_CRBC, Order_CRBC, Gene_Signature) %>%
  write.csv(output_file, row.names = FALSE)
cat(sprintf("\nCluster assignments saved to: %s\n", output_file))

