#!/usr/bin/env Rscript

# 0) Ensure required packages are available
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

required_pkgs <- c("argparse", "data.table", "dplyr", "tidyr", "ggplot2", "ggVennDiagram", "UpSetR")
invisible(lapply(required_pkgs, install_if_missing))

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggVennDiagram)
  library(UpSetR)
})

# 1) Parse arguments
parser <- ArgumentParser(
  description = "Create Venn diagrams and combination statistics for Chromosome/Plasmid/Virus by crop group."
)
parser$add_argument("-i", "--input", required = TRUE, help = "Input CSV file path")
parser$add_argument("-o", "--output", required = TRUE, help = "Output directory path")
parser$add_argument("--crop_col", default = "Host", help = "Crop grouping column name")
parser$add_argument("--contig_col", default = "Contig_Type3", help = "Contig type column name")
parser$add_argument(
  "--analysis_type",
  default = "defense",
  help = "Feature preset: defense, antidefense, amr"
)
parser$add_argument(
  "--feature_col",
  default = "",
  help = "Optional custom feature column; overrides --analysis_type"
)
parser$add_argument(
  "--crops",
  default = "",
  help = "Optional comma-separated crop names; if empty, use all crops"
)
parser$add_argument(
  "--four_crops",
  default = "",
  help = "Optional comma-separated 4 crop names for contig-based 4-set Venn; if empty, auto-select 4 crops"
)
args <- parser$parse_args()

# Resolve feature column by preset or custom override
analysis_type <- tolower(trimws(args$analysis_type))
preset_feature_map <- c(
  defense = "Defense_Subtype",
  antidefense = "AntiDS_Type",
  amr = "AMR_Type"
)

if (!(analysis_type %in% names(preset_feature_map))) {
  stop("Invalid --analysis_type. Allowed values: defense, antidefense, amr")
}

selected_feature_col <- if (trimws(args$feature_col) != "") {
  trimws(args$feature_col)
} else {
  unname(preset_feature_map[[analysis_type]])
}

# 2) Create output directory
if (!dir.exists(args$output)) {
  dir.create(args$output, recursive = TRUE)
}

# 3) Read input
dt <- fread(args$input, header = TRUE, stringsAsFactors = FALSE)

# 4) Validate required columns
required_cols <- c(args$crop_col, args$contig_col, selected_feature_col)
missing_cols <- setdiff(required_cols, colnames(dt))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}
cat(sprintf("Analysis type: %s | Feature column: %s\n", analysis_type, selected_feature_col))

# 5) Normalize contig groups to Chromosome/Plasmid/Virus
normalize_contig_group <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(
    grepl("plasmid", x),
    "Plasmid",
    ifelse(
      grepl("virus|phage|provirus", x),
      "Virus",
      ifelse(grepl("chromosome", x), "Chromosome", NA_character_)
    )
  )
}

work_df <- dt %>%
  mutate(
    Crop = .data[[args$crop_col]],
    Contig_Group = normalize_contig_group(.data[[args$contig_col]]),
    Feature_Raw = .data[[selected_feature_col]]
  ) %>%
  filter(!is.na(Crop), Crop != "", !is.na(Contig_Group)) %>%
  mutate(Crop = trimws(as.character(Crop)))

# 6) Expand comma-separated features
feature_df <- work_df %>%
  mutate(Feature_Raw = as.character(Feature_Raw)) %>%
  filter(!is.na(Feature_Raw), Feature_Raw != "") %>%
  separate_rows(Feature_Raw, sep = ",") %>%
  mutate(Feature = trimws(Feature_Raw)) %>%
  filter(!is.na(Feature), Feature != "") %>%
  distinct(Crop, Contig_Group, Feature)

# 7) Filter crops if provided
all_crops <- sort(unique(feature_df$Crop))
if (args$crops != "") {
  selected_crops <- trimws(unlist(strsplit(args$crops, ",")))
  selected_crops <- selected_crops[selected_crops != ""]
  all_crops <- intersect(all_crops, selected_crops)
}

if (length(all_crops) == 0) {
  stop("No crop groups available after filtering.")
}

# 8) Helper to compute all non-empty set combinations
compute_combination_table <- function(set_list, group_label_col, group_label_value) {
  set_names <- names(set_list)
  combo_rows <- list()
  idx <- 1L

  for (k in seq_along(set_names)) {
    combos <- combn(set_names, k, simplify = FALSE)
    for (cmb in combos) {
      inter_vals <- Reduce(intersect, set_list[cmb])
      non_cmb <- setdiff(set_names, cmb)
      union_non_cmb <- if (length(non_cmb) == 0) {
        character(0)
      } else {
        unique(unlist(set_list[non_cmb], use.names = FALSE))
      }
      exclusive_vals <- setdiff(inter_vals, union_non_cmb)

      combo_rows[[idx]] <- data.frame(
        Group_Label = group_label_value,
        Group_Label_Column = group_label_col,
        Combination = paste(cmb, collapse = "&"),
        Set_Count = k,
        Intersection_Size = length(inter_vals),
        Exclusive_Size = length(exclusive_vals),
        Intersection_Features = paste(sort(inter_vals), collapse = ";"),
        Exclusive_Features = paste(sort(exclusive_vals), collapse = ";"),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  bind_rows(combo_rows)
}

# Helper to sanitize names for output files
safe_name <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}

# Helper to build set mapping for UpSet outputs
build_upset_mapping <- function(set_list) {
  set_labels <- names(set_list)
  set_cols <- make.unique(paste0("Set_", seq_along(set_labels), "_", safe_name(set_labels)))
  data.frame(Set_Column = set_cols, Set_Label = set_labels, stringsAsFactors = FALSE)
}

# Helper to build binary matrix for UpSetR
build_upset_binary_matrix <- function(set_list) {
  universe <- sort(unique(unlist(set_list, use.names = FALSE)))
  if (length(universe) == 0) {
    return(NULL)
  }

  upset_df <- data.frame(Feature = universe, stringsAsFactors = FALSE)
  for (set_name in names(set_list)) {
    upset_df[[set_name]] <- as.integer(universe %in% set_list[[set_name]])
  }
  upset_df
}

# 9) Analysis A: per-crop, compare Chromosome/Plasmid/Virus
venn_dir_a <- file.path(args$output, "venn_by_crop_contig3")
if (!dir.exists(venn_dir_a)) {
  dir.create(venn_dir_a, recursive = TRUE)
}

set_size_rows <- list()
unique_rows <- list()
combo_rows <- list()
row_i <- 1L

for (crop_name in all_crops) {
  crop_df <- feature_df %>% filter(Crop == crop_name)

  set_list <- list(
    Chromosome = sort(unique(crop_df$Feature[crop_df$Contig_Group == "Chromosome"])),
    Plasmid = sort(unique(crop_df$Feature[crop_df$Contig_Group == "Plasmid"])),
    Virus = sort(unique(crop_df$Feature[crop_df$Contig_Group == "Virus"]))
  )

  # 9.1) Save set sizes
  for (set_name in names(set_list)) {
    set_size_rows[[length(set_size_rows) + 1L]] <- data.frame(
      Crop = crop_name,
      Contig_Group = set_name,
      Feature_Count = length(set_list[[set_name]]),
      stringsAsFactors = FALSE
    )
  }

  # 9.2) Save unique features per set
  for (set_name in names(set_list)) {
    other_names <- setdiff(names(set_list), set_name)
    other_union <- unique(unlist(set_list[other_names], use.names = FALSE))
    uniq_vals <- setdiff(set_list[[set_name]], other_union)
    if (length(uniq_vals) > 0) {
      unique_rows[[length(unique_rows) + 1L]] <- data.frame(
        Crop = crop_name,
        Contig_Group = set_name,
        Unique_Feature = uniq_vals,
        stringsAsFactors = FALSE
      )
    }
  }

  # 9.3) Save all combination statistics
  combo_rows[[row_i]] <- compute_combination_table(set_list, "Crop", crop_name)
  row_i <- row_i + 1L

  # 9.4) Draw Venn diagram
  venn_plot <- ggVennDiagram(
    x = set_list,
    label = "count",
    label_size = 5.2,
    edge_size = 0.9,
    set_size = 5.5
  ) +
    scale_fill_gradient(low = "#dce6f2", high = "#3f4f8a") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    ) +
    ggtitle(sprintf("%s: Chromosome vs Plasmid vs Virus", crop_name))

  safe_crop <- safe_name(crop_name)
  venn_pdf <- file.path(venn_dir_a, sprintf("venn_%s.pdf", safe_crop))
  ggsave(venn_pdf, venn_plot, width = 9, height = 7, device = "pdf")
}

# 10) Merge and save output tables for Analysis A
set_size_df <- bind_rows(set_size_rows) %>% arrange(Crop, Contig_Group)
combo_df <- bind_rows(combo_rows) %>%
  arrange(Group_Label, Set_Count, Combination) %>%
  rename(Crop = Group_Label)
unique_df <- if (length(unique_rows) > 0) {
  bind_rows(unique_rows) %>% arrange(Crop, Contig_Group, Unique_Feature)
} else {
  data.frame(
    Crop = character(),
    Contig_Group = character(),
    Unique_Feature = character(),
    stringsAsFactors = FALSE
  )
}

fwrite(set_size_df, file.path(args$output, "venn_set_sizes_by_crop_contig3.csv"))
fwrite(combo_df, file.path(args$output, "venn_all_combinations_by_crop_contig3.csv"))
fwrite(unique_df, file.path(args$output, "venn_unique_features_by_crop_contig3.csv"))

# 11) Analysis B: per-contig-group, compare 4 crops
venn_dir_b <- file.path(args$output, "venn_by_contig_crop4")
if (!dir.exists(venn_dir_b)) {
  dir.create(venn_dir_b, recursive = TRUE)
}
upset_dir_b <- file.path(args$output, "upset_by_contig_crop4")
if (!dir.exists(upset_dir_b)) {
  dir.create(upset_dir_b, recursive = TRUE)
}

crop_feature_counts <- feature_df %>%
  group_by(Crop) %>%
  summarise(Feature_Count = n_distinct(Feature), .groups = "drop") %>%
  arrange(desc(Feature_Count), Crop)

if (args$four_crops != "") {
  selected_4_crops <- trimws(unlist(strsplit(args$four_crops, ",")))
  selected_4_crops <- selected_4_crops[selected_4_crops != ""]
} else {
  selected_4_crops <- crop_feature_counts$Crop[1:min(4, nrow(crop_feature_counts))]
}

selected_4_crops <- unique(selected_4_crops)
selected_4_crops <- selected_4_crops[selected_4_crops %in% unique(feature_df$Crop)]

if (length(selected_4_crops) != 4) {
  stop("Exactly 4 crops are required for Analysis B. Please provide --four_crops with 4 valid crop names.")
}

contig_levels <- c("Chromosome", "Plasmid", "Virus")
set_size_rows_b <- list()
combo_rows_b <- list()
unique_rows_b <- list()
row_j <- 1L

for (contig_name in contig_levels) {
  contig_df <- feature_df %>%
    filter(Contig_Group == contig_name, Crop %in% selected_4_crops)

  set_list_b <- lapply(selected_4_crops, function(x) {
    sort(unique(contig_df$Feature[contig_df$Crop == x]))
  })
  names(set_list_b) <- selected_4_crops

  for (crop_name in names(set_list_b)) {
    set_size_rows_b[[length(set_size_rows_b) + 1L]] <- data.frame(
      Contig_Group = contig_name,
      Crop = crop_name,
      Feature_Count = length(set_list_b[[crop_name]]),
      stringsAsFactors = FALSE
    )
  }

  # Save unique features per crop under the same contig group
  for (crop_name in names(set_list_b)) {
    other_crops <- setdiff(names(set_list_b), crop_name)
    other_union <- unique(unlist(set_list_b[other_crops], use.names = FALSE))
    uniq_vals <- setdiff(set_list_b[[crop_name]], other_union)
    if (length(uniq_vals) > 0) {
      unique_rows_b[[length(unique_rows_b) + 1L]] <- data.frame(
        Contig_Group = contig_name,
        Crop = crop_name,
        Unique_Feature = uniq_vals,
        stringsAsFactors = FALSE
      )
    }
  }

  combo_rows_b[[row_j]] <- compute_combination_table(set_list_b, "Contig_Group", contig_name)
  row_j <- row_j + 1L

  venn_plot_b <- ggVennDiagram(
    x = set_list_b,
    label = "count",
    label_size = 4.8,
    edge_size = 0.9,
    set_size = 5.1
  ) +
    scale_fill_gradient(low = "#dce6f2", high = "#3f4f8a") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    ) +
    ggtitle(sprintf("%s: 4-Crop Venn", contig_name))

  venn_pdf_b <- file.path(venn_dir_b, sprintf("venn_%s_4crops.pdf", safe_name(contig_name)))
  ggsave(venn_pdf_b, venn_plot_b, width = 9, height = 7, device = "pdf")

  upset_matrix_b <- build_upset_binary_matrix(set_list_b)
  if (!is.null(upset_matrix_b)) {
    upset_pdf_b <- file.path(upset_dir_b, sprintf("upset_%s_4crops.pdf", safe_name(contig_name)))
    pdf(upset_pdf_b, width = 11, height = 7)
    upset_obj <- UpSetR::upset(
      upset_matrix_b,
      nsets = length(set_list_b),
      sets = names(set_list_b),
      order.by = "freq",
      keep.order = TRUE,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size"
    )
    print(upset_obj)
    dev.off()

    matrix_csv_b <- file.path(upset_dir_b, sprintf("upset_matrix_%s_4crops.csv", safe_name(contig_name)))
    fwrite(upset_matrix_b, matrix_csv_b)
  }

  mapping_csv_b <- file.path(upset_dir_b, sprintf("upset_set_mapping_%s_4crops.csv", safe_name(contig_name)))
  fwrite(build_upset_mapping(set_list_b), mapping_csv_b)
}

set_size_df_b <- bind_rows(set_size_rows_b) %>% arrange(Contig_Group, Crop)
combo_df_b <- bind_rows(combo_rows_b) %>%
  arrange(Group_Label, Set_Count, Combination) %>%
  rename(Contig_Group = Group_Label)
unique_df_b <- if (length(unique_rows_b) > 0) {
  bind_rows(unique_rows_b) %>% arrange(Contig_Group, Crop, Unique_Feature)
} else {
  data.frame(
    Contig_Group = character(),
    Crop = character(),
    Unique_Feature = character(),
    stringsAsFactors = FALSE
  )
}

fwrite(set_size_df_b, file.path(args$output, "venn_set_sizes_by_contig_crop4.csv"))
fwrite(combo_df_b, file.path(args$output, "venn_all_combinations_by_contig_crop4.csv"))
fwrite(unique_df_b, file.path(args$output, "venn_unique_features_by_contig_crop4.csv"))

# 12) Additional tidy outputs for exact combination statistics
exact_counts_a <- combo_df %>%
  transmute(
    Crop,
    Combination,
    Set_Count,
    Exact_Count = Exclusive_Size,
    Intersection_Size
  ) %>%
  arrange(Crop, Set_Count, Combination)
exact_features_a <- combo_df %>%
  select(Crop, Combination, Set_Count, Exclusive_Features) %>%
  filter(!is.na(Exclusive_Features), Exclusive_Features != "") %>%
  separate_rows(Exclusive_Features, sep = ";") %>%
  rename(Feature = Exclusive_Features) %>%
  arrange(Crop, Set_Count, Combination, Feature)

exact_counts_b <- combo_df_b %>%
  transmute(
    Contig_Group,
    Combination,
    Set_Count,
    Exact_Count = Exclusive_Size,
    Intersection_Size
  ) %>%
  arrange(Contig_Group, Set_Count, Combination)
exact_features_b <- combo_df_b %>%
  select(Contig_Group, Combination, Set_Count, Exclusive_Features) %>%
  filter(!is.na(Exclusive_Features), Exclusive_Features != "") %>%
  separate_rows(Exclusive_Features, sep = ";") %>%
  rename(Feature = Exclusive_Features) %>%
  arrange(Contig_Group, Set_Count, Combination, Feature)

fwrite(exact_counts_a, file.path(args$output, "venn_exact_combination_counts_by_crop_contig3.csv"))
fwrite(exact_features_a, file.path(args$output, "venn_exact_combination_features_by_crop_contig3.csv"))
fwrite(exact_counts_b, file.path(args$output, "venn_exact_combination_counts_by_contig_crop4.csv"))
fwrite(exact_features_b, file.path(args$output, "venn_exact_combination_features_by_contig_crop4.csv"))

# 13) Analysis C: per-crop, UpSet for Chromosome/Plasmid/Virus
upset_dir_c <- file.path(args$output, "upset_by_crop_contig3")
if (!dir.exists(upset_dir_c)) {
  dir.create(upset_dir_c, recursive = TRUE)
}

for (crop_name in all_crops) {
  crop_df <- feature_df %>% filter(Crop == crop_name)

  set_list_c <- list(
    Chromosome = sort(unique(crop_df$Feature[crop_df$Contig_Group == "Chromosome"])),
    Plasmid = sort(unique(crop_df$Feature[crop_df$Contig_Group == "Plasmid"])),
    Virus = sort(unique(crop_df$Feature[crop_df$Contig_Group == "Virus"]))
  )

  upset_matrix_c <- build_upset_binary_matrix(set_list_c)
  if (!is.null(upset_matrix_c)) {
    upset_pdf_c <- file.path(upset_dir_c, sprintf("upset_%s_contig3.pdf", safe_name(crop_name)))
    pdf(upset_pdf_c, width = 11, height = 7)
    upset_obj_c <- UpSetR::upset(
      upset_matrix_c,
      nsets = length(set_list_c),
      sets = names(set_list_c),
      order.by = "freq",
      keep.order = TRUE,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size"
    )
    print(upset_obj_c)
    dev.off()

    matrix_csv_c <- file.path(upset_dir_c, sprintf("upset_matrix_%s_contig3.csv", safe_name(crop_name)))
    fwrite(upset_matrix_c, matrix_csv_c)
  }

  mapping_csv_c <- file.path(upset_dir_c, sprintf("upset_set_mapping_%s_contig3.csv", safe_name(crop_name)))
  fwrite(build_upset_mapping(set_list_c), mapping_csv_c)
}

# Reuse Analysis A unique features for per-crop UpSet view
fwrite(unique_df, file.path(args$output, "upset_unique_features_by_crop_contig3.csv"))

# Reuse Analysis B unique features for per-contig UpSet view
fwrite(unique_df_b, file.path(args$output, "upset_unique_features_by_contig_crop4.csv"))

cat("Analysis completed.\n")
cat(sprintf("Output directory: %s\n", args$output))
cat("Generated files:\n")
cat("  - venn_set_sizes_by_crop_contig3.csv\n")
cat("  - venn_all_combinations_by_crop_contig3.csv\n")
cat("  - venn_unique_features_by_crop_contig3.csv\n")
cat("  - venn_by_crop_contig3/*.pdf\n")
cat("  - venn_set_sizes_by_contig_crop4.csv\n")
cat("  - venn_all_combinations_by_contig_crop4.csv\n")
cat("  - venn_unique_features_by_contig_crop4.csv\n")
cat("  - venn_by_contig_crop4/*.pdf\n")
cat("  - venn_exact_combination_counts_by_crop_contig3.csv\n")
cat("  - venn_exact_combination_features_by_crop_contig3.csv\n")
cat("  - venn_exact_combination_counts_by_contig_crop4.csv\n")
cat("  - venn_exact_combination_features_by_contig_crop4.csv\n")
cat("  - upset_by_contig_crop4/*.pdf\n")
cat("  - upset_by_contig_crop4/upset_set_mapping_*.csv\n")
cat("  - upset_by_contig_crop4/upset_matrix_*.csv\n")
cat("  - upset_unique_features_by_contig_crop4.csv\n")
cat("  - upset_by_crop_contig3/*.pdf\n")
cat("  - upset_by_crop_contig3/upset_set_mapping_*.csv\n")
cat("  - upset_by_crop_contig3/upset_matrix_*.csv\n")
cat("  - upset_unique_features_by_crop_contig3.csv\n")
