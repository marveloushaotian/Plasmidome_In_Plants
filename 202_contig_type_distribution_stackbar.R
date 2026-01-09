#!/usr/bin/env Rscript

# =============================================================================
# Contig Type Distribution Stacked Bar Charts
# Description: Generate stacked bar charts showing distribution of gene subtypes
#              (Defense_Subtype, AntiDS_Type, AMR_Type) by Host and Contig_Type
#              Creates both percentage and count visualizations
# Usage: Rscript 202_contig_type_distribution_stackbar.R -i <input.csv> -o <output_dir> [-n <top_n>]
#
# Arguments:
#   -i: Input CSV file path (required)
#   -o: Output base directory for results (required)
#   -n: Number of top subtypes to show (others merged into "Others") (default: 20)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(argparse)
})

# Parse command line arguments
parser <- ArgumentParser(description = "Generate stacked bar charts for gene subtype distribution")
parser$add_argument("-i", "--input", required = TRUE,
                    help = "Input CSV file path (e.g., Contig_Sample_Mapping_Final.csv)")
parser$add_argument("-o", "--output", required = TRUE,
                    help = "Output base directory for results")
parser$add_argument("-n", "--top-n", type = "integer", default = 20,
                    help = "Number of top subtypes to show, others merged to 'Others' (default: 20)")

args <- parser$parse_args()

# Fixed 21 colors for SubType
subtype_colors <- c(
  "#434d91","#6566aa","#9b7baa","#c6a4c5","#c6f0ec","#8fced1","#53a4a6","#d0cab7","#c0dbe6",
  "#509d95","#75b989","#92ca77","#d6ecc1","#e7ee9f","#f7ded5","#faaf7f","#f07e40","#dc5772",
  "#ebc1d1","#f9e7e7","#d1d9e2"
)

# Function to ensure directory exists
ensure_dir <- function(path) {
  if(!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", path))
  }
}

# Output subdirectories
output_dirs <- list(
  pct = file.path(args$output, "percentage"),
  count = file.path(args$output, "count")
)

# Create output directories
cat(sprintf("Output base directory: %s\n", args$output))
lapply(output_dirs, ensure_dir)

# Read data
cat(sprintf("Reading data from: %s\n", args$input))
data <- read.csv(args$input, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# 合并Virus和Provirus到Virus
data_processed <- data %>%
  mutate(
    Contig_Type = ifelse(Contig_Type %in% c("Virus", "Provirus"), "Virus", Contig_Type)
  ) %>%
  filter(Host != "" & Contig_Type != "")

# 要处理的三个子类型（反映需求：AntiDS_Subtype部分改为使用AntiDS_Type）
subtype_list <- list(
  Defense_Subtype = "Defense_Subtype",
  AntiDS_Subtype = "AntiDS_Type",
  AMR_Type = "AMR_Type"
)

# 定义Contig_Type和Host顺序与配色
contig_colors <- c(
  "Chromosome" = "#434d91",
  "Plasmid" = "#75b989",
  "Virus" = "#f07e40"
)
contig_order <- c("Chromosome", "Plasmid", "Virus")

host_order <- sort(unique(data_processed$Host))

# 需要去除的defense系统名称
defense_exclude_list <- c("VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA")

for (subtype in names(subtype_list)) {
  subtype_col <- subtype_list[[subtype]]
  # 判断该列是否存在于数据
  if (!(subtype_col %in% colnames(data_processed))) {
    message(paste0("Warning: Column ", subtype_col, " not found in the data, skipping..."))
    next
  }

  # 只保留必要的列
  df_sub <- data_processed %>%
    select(Host, Contig_Type, !!sym(subtype_col))

  # 把NA替换为""
  df_sub[[subtype_col]][is.na(df_sub[[subtype_col]])] <- ""

  # 分割并展开
  df_expanded <- df_sub %>%
    mutate(!!subtype_col := strsplit(as.character(!!sym(subtype_col)), ",")) %>%
    unnest(!!sym(subtype_col))

  # 去空白符和空行
  df_expanded[[subtype_col]] <- str_trim(df_expanded[[subtype_col]])
  df_expanded <- df_expanded %>%
    filter(!!sym(subtype_col) != "")

  # 关键修改：如果是Defense_Subtype，则去掉名称中包含"other"(大小写不敏感)以及明确排除的防御系统
  if (subtype == "Defense_Subtype") {
    df_expanded <- df_expanded %>%
      filter(
        !grepl("other", !!sym(subtype_col), ignore.case = TRUE),
        ! (!!sym(subtype_col) %in% defense_exclude_list)
      )
  }

  # 统计每种Subtype总体计数（全局）
  global_subtype_count <- df_expanded %>%
    group_by(Subtype = !!sym(subtype_col)) %>%
    summarise(Global_Count = n()) %>%
    arrange(desc(Global_Count))

  # Get top N subtypes based on global count
  top_n_subtypes <- global_subtype_count$Subtype[1:min(args$top_n, nrow(global_subtype_count))]
  all_subtypes <- global_subtype_count$Subtype
  cat(sprintf("  Found %d unique subtypes, selecting top %d\n", 
              length(all_subtypes), length(top_n_subtypes)))
  
  # Set final levels order: top N + Others (max N+1 categories)
  # Colors: highest abundance at bottom, Others at top
  final_subtype_order <- rev(c(top_n_subtypes, "Others"))

  # Mark "Others" for subtypes not in top N
  df_expanded <- df_expanded %>%
    mutate(
      Subtype = ifelse(!!sym(subtype_col) %in% top_n_subtypes, !!sym(subtype_col), "Others")
    )

  # 统计计数(已经合并Others)
  count_subtype <- df_expanded %>%
    group_by(Host, Contig_Type, Subtype) %>%
    summarise(Count = n(), .groups = 'drop')

  # 计算百分比 (以每个Host/Contig_Type为100%)
  percentage_subtype <- count_subtype %>%
    group_by(Host, Contig_Type) %>%
    mutate(
      Total = sum(Count),
      Percentage = (Count / Total) * 100
    ) %>%
    ungroup() %>%
    select(-Total) %>%
    # Ensure percentages sum to exactly 100% per group by adjusting rounding errors
    group_by(Host, Contig_Type) %>%
    mutate(
      Percentage = ifelse(Percentage > 100, 100, Percentage),  # Cap at 100
      Sum_Pct = sum(Percentage),
      Percentage = ifelse(Sum_Pct > 100, Percentage * (100 / Sum_Pct), Percentage)
    ) %>%
    select(-Sum_Pct) %>%
    ungroup()

  # 水平与配色
  plot_data_pct <- percentage_subtype %>%
    mutate(
      Host = factor(Host, levels = host_order),
      Contig_Type = factor(Contig_Type, levels = rev(contig_order)),
      Subtype = factor(Subtype, levels = final_subtype_order)
    )
  plot_data_count <- count_subtype %>%
    mutate(
      Host = factor(Host, levels = host_order),
      Contig_Type = factor(Contig_Type, levels = rev(contig_order)),
      Subtype = factor(Subtype, levels = final_subtype_order)
    )

  # 固定21色（与order一一对应）, 同样反转颜色顺序
  used_colors <- subtype_colors[seq_len(length(final_subtype_order))]
  names(used_colors) <- rev(final_subtype_order)  # Important for matching: lowest is lowest, Others is top

  # --- 保存百分比和计数的plot与表到不同子目录 ---

  for (ct in contig_order) {
    plot_pct <- plot_data_pct %>%
      filter(Contig_Type == ct)
    if(nrow(plot_pct) > 0) {
      # 图形中堆叠顺序仍然用final_subtype_order（和左图条形堆叠一致，最下颜色是丰度最高）
      # 图例顺序反转：使用 override.aes = list(order = ...) 来手动调整图例顺序
      # ggplot2 没有直接的图例顺序参数，最佳方法是更改factor顺序，仅用于图例
      # 我们为legend单独做一个反向的factor
      plot_pct$Subtype_legend <- factor(plot_pct$Subtype, levels = rev(final_subtype_order))

      p_pct <- ggplot(plot_pct, aes(x=Host, y=Percentage, fill=Subtype)) +
        geom_bar(stat="identity", width=0.7) +
        scale_fill_manual(
          values=used_colors,
          drop=FALSE,
          # 指定breaks为反向顺序，只影响图例顺序，不影响左边显示顺序
          breaks = rev(final_subtype_order)
        ) +
        theme_bw() +
        theme(
          axis.text.x = element_text(hjust=0.5, size=13),
          axis.text.y = element_text(size=13),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.position="right"
        ) +
        guides(
          fill = guide_legend(ncol = 1, reverse = TRUE) # 反转图例顺序
        ) +
        labs(
          y="Percentage (%)",
          fill=paste0(subtype, " (", ct, ")")
        ) +
        scale_y_continuous(expand = c(0,0), limits = c(0,100), oob = scales::squish)
      ggsave(
        file.path(output_dirs[["pct"]], sprintf("%s_distribution_percentage_byhost_%s.pdf", tolower(subtype), tolower(ct))),
        p_pct, width=10, height=7
      )
    }
  }

  for (ct in contig_order) {
    plot_count <- plot_data_count %>%
      filter(Contig_Type == ct)
    if(nrow(plot_count) > 0) {
      # 同样调整图例顺序
      plot_count$Subtype_legend <- factor(plot_count$Subtype, levels = rev(final_subtype_order))
      p_count <- ggplot(plot_count, aes(x=Host, y=Count, fill=Subtype)) +
        geom_bar(stat="identity", width=0.7) +
        scale_fill_manual(
          values=used_colors,
          drop=FALSE,
          breaks = rev(final_subtype_order)
        ) +
        theme_bw() +
        theme(
          axis.text.x = element_text(hjust=0.5, size=13),
          axis.text.y = element_text(size=13),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.position="right"
        ) +
        guides(
          fill = guide_legend(ncol = 1, reverse = TRUE)
        ) +
        labs(
          y="Count",
          fill=paste0(subtype, " (", ct, ")")
        ) +
        scale_y_continuous(expand=c(0,0))
      ggsave(
        file.path(output_dirs[["count"]], sprintf("%s_distribution_count_byhost_%s.pdf", tolower(subtype), tolower(ct))),
        p_count, width=10, height=7
      )
    }
  }
  # 保存数据
  write.csv(
    count_subtype,
    file.path(output_dirs[["count"]], sprintf("%s_distribution_count.csv", tolower(subtype))),
    row.names=FALSE
  )
  write.csv(
    percentage_subtype,
    file.path(output_dirs[["pct"]], sprintf("%s_distribution_percentage.csv", tolower(subtype))),
    row.names=FALSE
  )
}

cat("\n========================================\n")
cat("Analysis completed successfully!\n")
cat("========================================\n")
cat(sprintf("\nOutput directory: %s\n", args$output))
cat(sprintf("Generated files:\n"))
for (subtype in names(subtype_list)) {
  cat(sprintf("  %s:\n", subtype))
  for (ct in contig_order) {
    cat(sprintf("    - percentage/%s_distribution_percentage_byhost_%s.pdf\n", tolower(subtype), tolower(ct)))
    cat(sprintf("    - count/%s_distribution_count_byhost_%s.pdf\n", tolower(subtype), tolower(ct)))
  }
  cat(sprintf("    - count/%s_distribution_count.csv\n", tolower(subtype)))
  cat(sprintf("    - percentage/%s_distribution_percentage.csv\n", tolower(subtype)))
}
cat("\n")

