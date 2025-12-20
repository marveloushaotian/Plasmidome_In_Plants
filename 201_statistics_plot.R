# Statistics and plots for Defense_Subtype, AntiDS_Type, AMR_Type
# Top 100 subtypes with unified x-axis across facets (all bars in gray), with defense filter for Defense_Subtype and AntiDS_Type, 
# and also exclude 'Other' for Defense_Subtype and AntiDS_Type

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# Read CSV file
data <- read.csv('Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final.csv', stringsAsFactors = FALSE)

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

# Function to process and plot top 100 subtypes, all bars in gray, with optional defense filtering
process_and_plot <- function(df, col_name, output_prefix, defense_filter = FALSE) {
  # Count subtypes
  subtype_counts <- count_subtypes(df, col_name, defense_filter)
  
  # Get top 100 subtypes based on total count
  top_subtypes <- subtype_counts %>%
    group_by(subtype) %>%
    summarise(total_count = sum(count), .groups = 'drop') %>%
    arrange(desc(total_count)) %>%
    slice_head(n = 100) %>%
    pull(subtype)
  
  # Keep only top 100 subtypes
  counts_top100 <- subtype_counts %>%
    filter(subtype %in% top_subtypes)
  
  # Fill in missing combinations with zeros
  all_combos <- expand.grid(
    Host = unique(counts_top100$Host),
    Contig_Type = unique(counts_top100$Contig_Type),
    subtype = top_subtypes,
    stringsAsFactors = FALSE
  )
  
  counts_top100_full <- all_combos %>%
    left_join(counts_top100, by = c("Host", "Contig_Type", "subtype")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  # Set subtype as factor with descending order
  subtype_levels <- counts_top100 %>%
    group_by(subtype) %>%
    summarise(total = sum(count), .groups = 'drop') %>%
    arrange(desc(total)) %>%
    pull(subtype)
  
  counts_top100_full$subtype <- factor(counts_top100_full$subtype, levels = subtype_levels)
  
  # Save table
  write_csv(counts_top100_full, 
            paste0('Result/NCBI_4395_Batch/', 
                   output_prefix, '_Top100_Counts_ByGroup.csv'))
  
  # Plot (all gray)
  p <- ggplot(counts_top100_full, aes(x = subtype, y = count)) +
    geom_bar(stat = 'identity', position = 'dodge', fill = "gray60") +
    facet_wrap(~Host + Contig_Type, scales = "fixed", ncol = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste0('Top 100 ', col_name, ' Counts by Host and Contig Type'),
         x = col_name, y = 'Count')
  
  ggsave(
    paste0('Result/NCBI_4395_Batch/', 
           output_prefix, '_Top100_Counts_ByGroup.pdf'),
    plot = p,
    width = 20, height = 30,
    limitsize = FALSE
  )
  
  return(counts_top100_full)
}

# Process Defense_Subtype (with defense filter, also remove "Other")
cat("Processing Defense_Subtype...\n")
defense_results <- process_and_plot(chl, 'Defense_Subtype', 'Defense_Subtype', defense_filter = TRUE)

# Process AntiDS_Type (with defense filter, also remove "Other")
cat("Processing AntiDS_Type...\n")
antids_results <- process_and_plot(chl, 'AntiDS_Type', 'AntiDS_Type', defense_filter = TRUE)

# Process AMR_Type (no filter)
cat("Processing AMR_Type...\n")
amr_results <- process_and_plot(chl, 'AMR_Type', 'AMR_Type', defense_filter = FALSE)

cat("All processing completed!\n")

