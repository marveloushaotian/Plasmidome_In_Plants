#!/usr/bin/env Rscript

# =============================================================================
# Analyze AMR-Taxonomy Association
# Description: Investigate which bacterial groups preferentially carry specific AMR genes
# Methods: Chi-square test, enrichment analysis, heatmap visualization
# Usage: Rscript 215_analyze_amr_taxonomy_association.R -c <coords.csv> -m <matrix.csv> -d <workdir> -o <output_prefix>
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Analyze AMR-Taxonomy Association")
parser$add_argument("-c", "--coords", required = TRUE, 
                    help = "PCoA coordinates CSV file")
parser$add_argument("-m", "--matrix", required = TRUE,
                    help = "AMR matrix CSV file")
parser$add_argument("-d", "--workdir", default = ".",
                    help = "Working directory (default: current directory)")
parser$add_argument("-o", "--output-prefix", default = "AMR",
                    help = "Output file prefix (default: AMR)")

args <- parser$parse_args()

setwd(args$workdir)

cat("=== AMR-Taxonomy Association Analysis ===\n\n")

# Read data
coords <- read.csv(args$coords, stringsAsFactors = FALSE)
amr_matrix <- read.csv(args$matrix, stringsAsFactors = FALSE)

# Merge data
amr_data <- amr_matrix
rownames(amr_data) <- amr_data$Sample_ID
amr_data$Sample_ID <- NULL

common_samples <- intersect(coords$Sample_ID, rownames(amr_data))
coords_m <- coords %>% filter(Sample_ID %in% common_samples)
amr_m <- amr_data[coords_m$Sample_ID, ]

# Convert to binary (presence/absence)
amr_binary <- (amr_m > 0) * 1

# Filter AMR genes present in at least 2% of samples
amr_prevalence <- colMeans(amr_binary) * 100
keep_amr <- names(amr_prevalence[amr_prevalence >= 2])
amr_binary <- amr_binary[, keep_amr]

cat(sprintf("Samples: %d\n", nrow(coords_m)))
cat(sprintf("AMR genes (>=2%% samples): %d\n\n", length(keep_amr)))

# Method 1: Chi-square test for each AMR-Class combination
cat("=== Method 1: Chi-square Association Test ===\n\n")

classes <- unique(coords_m$Class_CRBC)
classes <- classes[!is.na(classes) & classes != ""]

results_chisq <- data.frame()

for (amr in colnames(amr_binary)) {
  for (cls in classes) {
    is_class <- coords_m$Class_CRBC == cls
    has_amr <- amr_binary[, amr] == 1
    
    n_class_amr <- sum(is_class & has_amr, na.rm = TRUE)
    n_class_no_amr <- sum(is_class & !has_amr, na.rm = TRUE)
    n_other_amr <- sum(!is_class & has_amr, na.rm = TRUE)
    n_other_no_amr <- sum(!is_class & !has_amr, na.rm = TRUE)
    
    cont_table <- matrix(c(n_class_amr, n_class_no_amr, n_other_amr, n_other_no_amr), 
                         nrow = 2, byrow = TRUE)
    
    if (min(cont_table) >= 5) {
      test <- chisq.test(cont_table, correct = FALSE)
      pval <- test$p.value
    } else {
      test <- fisher.test(cont_table)
      pval <- test$p.value
    }
    
    if (n_class_no_amr > 0 & n_other_amr > 0 & n_other_no_amr > 0) {
      odds_ratio <- (n_class_amr / n_class_no_amr) / (n_other_amr / n_other_no_amr)
    } else {
      odds_ratio <- NA
    }
    
    prev_class <- n_class_amr / (n_class_amr + n_class_no_amr) * 100
    prev_overall <- sum(has_amr, na.rm = TRUE) / length(has_amr) * 100
    
    results_chisq <- rbind(results_chisq, data.frame(
      AMR = amr,
      Class = cls,
      N_Class_with_AMR = n_class_amr,
      N_Class_total = n_class_amr + n_class_no_amr,
      Prevalence_Class = prev_class,
      Prevalence_Overall = prev_overall,
      Odds_Ratio = odds_ratio,
      P_value = pval
    ))
  }
}

# Adjust p-values (FDR)
results_chisq$P_adjusted <- p.adjust(results_chisq$P_value, method = "BH")
results_chisq$Log2_OR <- log2(results_chisq$Odds_Ratio)

# Filter significant associations
sig_assoc <- results_chisq %>%
  filter(P_adjusted < 0.05 & !is.na(Odds_Ratio)) %>%
  arrange(P_adjusted)

cat(sprintf("Significant associations (FDR < 0.05): %d\n\n", nrow(sig_assoc)))

# Top enrichments
cat("=== Top Enriched AMR-Class Associations (Top 20) ===\n")
top_enriched <- sig_assoc %>%
  filter(Odds_Ratio > 1) %>%
  arrange(desc(Odds_Ratio)) %>%
  head(20)

print(top_enriched %>% select(Class, AMR, Prevalence_Class, Prevalence_Overall, Odds_Ratio, P_adjusted))

# Top depletions
cat("\n=== Top Depleted AMR-Class Associations (Top 10) ===\n")
top_depleted <- sig_assoc %>%
  filter(Odds_Ratio < 1) %>%
  arrange(Odds_Ratio) %>%
  head(10)

print(top_depleted %>% select(Class, AMR, Prevalence_Class, Prevalence_Overall, Odds_Ratio, P_adjusted))

# Save full results
output_file <- sprintf("%s_Class_Association_ChiSquare.csv", args$output_prefix)
write.csv(results_chisq, output_file, row.names = FALSE)
cat(sprintf("\nFull results saved: %s\n", output_file))

# Method 2: Heatmap of AMR prevalence by Class
cat("\n=== Method 2: AMR Prevalence Heatmap ===\n")

amr_by_class <- coords_m %>%
  select(Sample_ID, Class_CRBC) %>%
  cbind(amr_binary) %>%
  group_by(Class_CRBC) %>%
  summarise(across(starts_with("AMR_"), ~mean(.) * 100), .groups = "drop")

prev_matrix <- as.matrix(amr_by_class[, -1])
rownames(prev_matrix) <- amr_by_class$Class_CRBC

amr_var <- apply(prev_matrix, 2, var)
top_var_amr <- names(sort(amr_var, decreasing = TRUE))[1:min(30, ncol(prev_matrix))]
prev_matrix_filt <- prev_matrix[, top_var_amr]

colnames(prev_matrix_filt) <- gsub("AMR_", "", colnames(prev_matrix_filt))

pdf(sprintf("%s_Class_Prevalence_Heatmap.pdf", args$output_prefix), width = 14, height = 8)
pheatmap(prev_matrix_filt,
         main = "AMR Gene Prevalence by Bacterial Class (%)",
         color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         fontsize_col = 8,
         angle_col = 45,
         border_color = NA,
         display_numbers = TRUE,
         number_format = "%.1f",
         fontsize_number = 6)
dev.off()
cat(sprintf("Heatmap saved: %s_Class_Prevalence_Heatmap.pdf\n", args$output_prefix))

# Method 3: Enrichment bubble plot
cat("\n=== Method 3: Enrichment Bubble Plot ===\n")

bubble_data <- sig_assoc %>%
  filter(Odds_Ratio > 1 & Prevalence_Class >= 5) %>%
  mutate(
    AMR_clean = gsub("AMR_", "", AMR),
    Significance = -log10(P_adjusted)
  )

if (nrow(bubble_data) > 0) {
  top_amr <- bubble_data %>%
    group_by(AMR_clean) %>%
    summarise(n_sig = n()) %>%
    arrange(desc(n_sig)) %>%
    head(15) %>%
    pull(AMR_clean)
  
  bubble_data_filt <- bubble_data %>%
    filter(AMR_clean %in% top_amr)
  
  p_bubble <- ggplot(bubble_data_filt, aes(x = AMR_clean, y = Class)) +
    geom_point(aes(size = Prevalence_Class, color = Log2_OR)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          name = "Log2(Odds Ratio)") +
    scale_size_continuous(range = c(2, 12), name = "Prevalence (%)") +
    labs(
      title = "AMR Gene Enrichment in Bacterial Classes",
      subtitle = "Only showing significant enrichments (FDR < 0.05, Prevalence >= 5%)",
      x = "AMR Gene",
      y = "Bacterial Class"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  ggsave(sprintf("%s_Class_Enrichment_Bubble.pdf", args$output_prefix), 
         p_bubble, width = 14, height = 8)
  cat(sprintf("Bubble plot saved: %s_Class_Enrichment_Bubble.pdf\n", args$output_prefix))
}

# Method 4: Summary statistics by Class
cat("\n=== Method 4: Summary Statistics by Class ===\n\n")

for (cls in sort(classes)) {
  cls_samples <- coords_m$Sample_ID[coords_m$Class_CRBC == cls]
  cls_amr <- amr_binary[cls_samples, , drop = FALSE]
  
  cls_prev <- colMeans(cls_amr) * 100
  top_amr <- sort(cls_prev[cls_prev > 0], decreasing = TRUE)[1:5]
  
  enriched <- sig_assoc %>%
    filter(Class == cls & Odds_Ratio > 2) %>%
    arrange(desc(Odds_Ratio)) %>%
    head(3)
  
  cat(sprintf("【%s】 (n=%d)\n", cls, length(cls_samples)))
  cat("  Most common AMR: ")
  cat(paste(sprintf("%s(%.1f%%)", gsub("AMR_", "", names(top_amr)), top_amr), collapse = ", "))
  cat("\n")
  
  if (nrow(enriched) > 0) {
    cat("  Specific enrichment: ")
    cat(paste(sprintf("%s(OR=%.1f)", gsub("AMR_", "", enriched$AMR), enriched$Odds_Ratio), collapse = ", "))
    cat("\n")
  }
  cat("\n")
}

# Method 5: Order-level analysis
cat("=== Method 5: Order-level Analysis ===\n\n")

orders <- unique(coords_m$Order_CRBC)
orders <- orders[!is.na(orders) & orders != ""]

order_counts <- table(coords_m$Order_CRBC)
major_orders <- names(order_counts[order_counts >= 30])

cat(sprintf("Major Orders (>=30 samples): %d\n\n", length(major_orders)))

results_order <- data.frame()

for (amr in colnames(amr_binary)) {
  for (ord in major_orders) {
    is_order <- coords_m$Order_CRBC == ord
    has_amr <- amr_binary[, amr] == 1
    
    n_order_amr <- sum(is_order & has_amr, na.rm = TRUE)
    n_order_no_amr <- sum(is_order & !has_amr, na.rm = TRUE)
    n_other_amr <- sum(!is_order & has_amr, na.rm = TRUE)
    n_other_no_amr <- sum(!is_order & !has_amr, na.rm = TRUE)
    
    cont_table <- matrix(c(n_order_amr, n_order_no_amr, n_other_amr, n_other_no_amr), 
                         nrow = 2, byrow = TRUE)
    
    if (min(cont_table) >= 5) {
      pval <- chisq.test(cont_table, correct = FALSE)$p.value
    } else {
      pval <- fisher.test(cont_table)$p.value
    }
    
    if (n_order_no_amr > 0 & n_other_amr > 0 & n_other_no_amr > 0) {
      odds_ratio <- (n_order_amr / n_order_no_amr) / (n_other_amr / n_other_no_amr)
    } else {
      odds_ratio <- NA
    }
    
    prev_order <- n_order_amr / (n_order_amr + n_order_no_amr) * 100
    
    results_order <- rbind(results_order, data.frame(
      AMR = amr,
      Order = ord,
      N_with_AMR = n_order_amr,
      N_total = n_order_amr + n_order_no_amr,
      Prevalence = prev_order,
      Odds_Ratio = odds_ratio,
      P_value = pval
    ))
  }
}

results_order$P_adjusted <- p.adjust(results_order$P_value, method = "BH")

top_order <- results_order %>%
  filter(P_adjusted < 0.01 & Odds_Ratio > 2 & Prevalence >= 10) %>%
  arrange(desc(Odds_Ratio)) %>%
  head(20)

cat("Top Order-level Enrichments (Top 20):\n")
print(top_order %>% select(Order, AMR, Prevalence, Odds_Ratio, P_adjusted))

write.csv(results_order, sprintf("%s_Order_Association.csv", args$output_prefix), row.names = FALSE)
cat(sprintf("\nOrder-level results saved: %s_Order_Association.csv\n", args$output_prefix))

cat("\n\n========================================\n")
cat("Analysis Complete! Output files:\n")
cat("========================================\n")
cat(sprintf("1. %s_Class_Association_ChiSquare.csv - Class-level association test\n", args$output_prefix))
cat(sprintf("2. %s_Class_Prevalence_Heatmap.pdf - AMR prevalence heatmap\n", args$output_prefix))
cat(sprintf("3. %s_Class_Enrichment_Bubble.pdf - Enrichment bubble plot\n", args$output_prefix))
cat(sprintf("4. %s_Order_Association.csv - Order-level association test\n", args$output_prefix))

