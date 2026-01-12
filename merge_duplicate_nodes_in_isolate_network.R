#!/usr/bin/env Rscript
# Merge duplicate rows by id column
# For label, AntiDS_Type, AMR_Type: concatenate non-empty values with semicolon
# For Defense_Num, AntiDS_Num, AMR_Num: sum the values

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript merge_duplicate_nodes_in_isolate_network.R input.tsv [output.tsv]")
}

input_file <- args[1]
if (length(args) >= 2) {
  output_file <- args[2]
} else {
  # Default output: remove _dup suffix if present, or add _merged suffix
  output_file <- sub("_dup\\.tsv$", ".tsv", input_file)
  if (output_file == input_file) {
    output_file <- sub("\\.tsv$", "_merged.tsv", input_file)
    if (output_file == input_file) {
      output_file <- paste0(input_file, "_merged.tsv")
    }
  }
}

cat(sprintf("Reading file: %s\n", input_file))
df <- read.delim(input_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Original rows: %d\n", nrow(df)))
cat(sprintf("Unique IDs: %d\n", length(unique(df$id))))

# Clean label column before merging
# Remove specific defense types: "VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA" and values ending with "other"
if ("label" %in% colnames(df)) {
  cat("Cleaning label column...\n")
  
  # Define values to remove
  values_to_remove <- c("VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA")
  
  # Function to clean a single label value
  clean_label <- function(label_str) {
    if (is.na(label_str) || label_str == "" || trimws(label_str) == "") {
      return("")
    }
    
    # Split by semicolon if multiple values
    values <- strsplit(label_str, ";")[[1]]
    values <- trimws(values)
    
    # Remove empty values
    values <- values[values != ""]
    
    # Remove specific values
    values <- values[!values %in% values_to_remove]
    
    # Remove values ending with "other" (case-insensitive)
    values <- values[!grepl("other$", values, ignore.case = TRUE)]
    
    # Return cleaned string
    if (length(values) == 0) {
      return("")
    } else {
      return(paste(values, collapse = ";"))
    }
  }
  
  # Apply cleaning to all label values
  df$label <- sapply(df$label, clean_label)
  
  cat("Label column cleaned.\n")
}

# Check required columns
required_cols <- c("id", "label", "AntiDS_Type", "AMR_Type", "Defense_Num", "AntiDS_Num", "AMR_Num")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s\n", paste(missing_cols, collapse = ", ")))
}

# Function to merge non-empty strings with semicolon
merge_strings <- function(x) {
  # Remove NA and empty strings
  non_empty <- x[!is.na(x) & x != "" & trimws(x) != ""]
  if (length(non_empty) == 0) {
    return("")
  }
  # Remove duplicates and merge
  unique_vals <- unique(non_empty)
  paste(unique_vals, collapse = ";")
}

# Function to sum numeric values
sum_numeric <- function(x) {
  # Convert to numeric, handling NA
  num_vals <- as.numeric(x)
  num_vals[is.na(num_vals)] <- 0
  sum(num_vals)
}

cat("Merging duplicate rows...\n")

# Group by id and aggregate
result_list <- list()

# Get all unique IDs
unique_ids <- unique(df$id)

# Process each ID
for (i in seq_along(unique_ids)) {
  id_val <- unique_ids[i]
  id_rows <- df[df$id == id_val, , drop = FALSE]
  
  # For string columns (label, AntiDS_Type, AMR_Type): merge with semicolon
  merged_row <- id_rows[1, , drop = FALSE]  # Start with first row
  
  # Merge label
  if ("label" %in% colnames(id_rows)) {
    merged_row$label <- merge_strings(id_rows$label)
  }
  
  # Merge AntiDS_Type
  if ("AntiDS_Type" %in% colnames(id_rows)) {
    merged_row$AntiDS_Type <- merge_strings(id_rows$AntiDS_Type)
  }
  
  # Merge AMR_Type
  if ("AMR_Type" %in% colnames(id_rows)) {
    merged_row$AMR_Type <- merge_strings(id_rows$AMR_Type)
  }
  
  # Sum numeric columns (Defense_Num, AntiDS_Num, AMR_Num)
  if ("Defense_Num" %in% colnames(id_rows)) {
    merged_row$Defense_Num <- sum_numeric(id_rows$Defense_Num)
  }
  
  if ("AntiDS_Num" %in% colnames(id_rows)) {
    merged_row$AntiDS_Num <- sum_numeric(id_rows$AntiDS_Num)
  }
  
  if ("AMR_Num" %in% colnames(id_rows)) {
    merged_row$AMR_Num <- sum_numeric(id_rows$AMR_Num)
  }
  
  result_list[[i]] <- merged_row
  
  if (i %% 1000 == 0) {
    cat(sprintf("Processed %d/%d unique IDs...\n", i, length(unique_ids)))
  }
}

# Combine results
result_df <- do.call(rbind, result_list)

cat(sprintf("Merged rows: %d\n", nrow(result_df)))
cat(sprintf("Reduction: %d rows removed\n", nrow(df) - nrow(result_df)))

# Write output
cat(sprintf("Writing output to: %s\n", output_file))
write.table(result_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Merge completed successfully!\n")
cat(sprintf("Output saved to: %s\n", output_file))
