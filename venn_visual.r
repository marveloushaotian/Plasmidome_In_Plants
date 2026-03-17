# Load required libraries
library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to create a list of unique Defense_Types for each group
create_defense_type_list <- function(data, group_columns) {
  # Group data by specified columns and summarize unique Defense_Types
  grouped_data <- data %>%
    group_by(across(all_of(group_columns))) %>%
    summarise(Defense_Types = list(unique(Defense_Subtype))) %>%
    ungroup()
  
  # Convert the grouped data into a named list
  defense_type_list <- setNames(grouped_data$Defense_Types, apply(grouped_data[group_columns], 1, paste, collapse = "_"))
  
  return(defense_type_list)
}

# Read the input CSV file
data <- read.csv("all_contig_combined_defense/r_all_contigs_info_more_defense_only.csv")

# Filter out rows where Contig_Classification is "Phage"
# data <- data %>% filter(Contig_Classification != "Phage")

# Specify the columns to group by
group_columns <- "Country" # Modify as needed

# Create the list of Defense_Types for each group
defense_type_list <- create_defense_type_list(data, group_columns)

# Create Venn diagram
venn_plot <- ggVennDiagram(
  label_alpha = 0,
  label_percent_digit = 2,  # Change to display percentages with two decimal places
  label_size = 9,
  defense_type_list,
  category.names = names(defense_type_list),
  label = "percent",
  edge_size = 0.5,
  set_size = 12
) +
  scale_fill_gradient(low = "#e8e7e9", high = "#6566aa") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "none")

# Save the plot as a PDF file
output_filename <- paste0("Results/statistic_venn/venn_", group_columns, ".pdf")
ggsave(output_filename, venn_plot, width = 10, height = 8, device = "pdf")

# Extract unique defense types for each group
unique_defenses <- list()
group_names <- names(defense_type_list)

# Find defense types that are unique to each group
for (i in 1:length(group_names)) {
  current_group <- group_names[i]
  current_defenses <- defense_type_list[[current_group]]
  
  # Find defenses that only appear in this group
  unique_to_current <- current_defenses
  for (j in 1:length(group_names)) {
    if (i != j) {
      unique_to_current <- setdiff(unique_to_current, defense_type_list[[group_names[j]]])
    }
  }
  
  unique_defenses[[current_group]] <- unique_to_current
}

# Create a data frame of unique defenses for each group
unique_defense_df <- data.frame(
  Group = character(),
  Unique_Defense = character(),
  stringsAsFactors = FALSE
)

for (group in names(unique_defenses)) {
  if (length(unique_defenses[[group]]) > 0) {
    temp_df <- data.frame(
      Group = rep(group, length(unique_defenses[[group]])),
      Unique_Defense = unique_defenses[[group]],
      stringsAsFactors = FALSE
    )
    unique_defense_df <- rbind(unique_defense_df, temp_df)
  }
}

# Save the unique defenses to a CSV file
unique_defense_output <- paste0("Results/statistic_venn/unique_defenses_", group_columns, ".csv")
write.csv(unique_defense_df, unique_defense_output, row.names = FALSE)

# Print summary of unique defenses
cat("Unique defense types for each group have been saved to:", unique_defense_output, "\n")
for (group in names(unique_defenses)) {
  cat("\nGroup:", group, "\n")
  cat("Number of unique defense types:", length(unique_defenses[[group]]), "\n")
  if (length(unique_defenses[[group]]) > 0) {
    cat("Unique defense types:", paste(unique_defenses[[group]], collapse = ", "), "\n")
  }
}
