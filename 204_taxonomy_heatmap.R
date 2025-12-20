# Generate heatmaps for gene distribution by taxonomy
# X轴: Genus (gtdb), Y轴: gene subtype（仅top50），Host为facet
# 只取每种基因值(top50)，超过部分不显示
# Genus (gtdb) 超过100个的也不显示
# 对Defense做指定过滤：去除特定defense系统以及包含"other"单词的defense（与stackbar脚本保持一致）

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)
library(stringr)

# 需要排除的defense系统
defense_exclude_list <- c("VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA")

# 读取数据
data <- read.csv('Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final.csv', 
                 stringsAsFactors = FALSE, check.names = FALSE)

# 合并Contig_Type
data <- data %>%
  mutate(Contig_Type = if_else(Contig_Type %in% c('Virus', 'Provirus'), 'Virus', Contig_Type))

# 基因计数函数，统计每个 Host+Contig_Type+Genus+gene subtype 的数量（支持defense过滤，排除名称含"other"的子类型）
count_by_genus <- function(df, gene_col, defense_filter = FALSE) {
  gene_df <- df %>%
    select(Host, Contig_Type, `Genus (gtdb)`, all_of(gene_col)) %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    mutate(gene = strsplit(.data[[gene_col]], ',')) %>%
    unnest(gene) %>%
    mutate(gene = trimws(gene)) %>%
    filter(gene != "", !is.na(`Genus (gtdb)`), `Genus (gtdb)` != "")
  if (defense_filter) {
    gene_df <- gene_df %>%
      filter(!(gene %in% defense_exclude_list)) %>%
      filter(!str_detect(tolower(gene), "other"))
  }
  gene_df %>%
    count(Host, Contig_Type, `Genus (gtdb)`, gene, name = 'count')
}

# 判断需要过滤的类型（defense 和 AntiDS 也要过滤含other的）
is_defense_like <- function(gene_col) {
  # 只有 Defense_Subtype 返回TRUE
  return(gene_col == "Defense_Subtype")
}

# 主绘图函数：每个Contig_Type分开，Host做facet，X为Genus，Y为gene（只显示top50 gene, Genus (gtdb)只显示最多100个）
create_heatmaps_by_contigtype <- function(df, gene_col, gene_name) {
  cat(paste0("Processing ", gene_name, "...\n"))
  
  # 基因计数（defense类做过滤，包括other；AntiDS不做过滤）
  gene_counts <- count_by_genus(df, gene_col, defense_filter = is_defense_like(gene_col))
  
  # 得到所有Contig_Type
  contig_types <- unique(gene_counts$Contig_Type)

  for (contig_type in contig_types) {
    cat(paste0("  Plotting for Contig_Type: ", contig_type, "\n"))
    
    # 筛选出当前Contig_Type的数据
    subdata <- gene_counts %>% filter(Contig_Type == contig_type)
    if (nrow(subdata) == 0) {
      cat(paste0("    No data for Contig_Type ", contig_type, ", skipping...\n"))
      next
    }

    # 计算该Contig_Type下，所有gene的总数，取top 50
    top_genes <- subdata %>%
      group_by(gene) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = 50) %>%
      pull(gene)

    # 只保留top 50 gene
    plot_data <- subdata %>%
      filter(gene %in% top_genes)

    # 只保留出现次数最多的前100个Genus (gtdb)
    genus_levels_full <- plot_data %>%
      group_by(`Genus (gtdb)`) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = 100) %>%
      pull(`Genus (gtdb)`)

    plot_data <- plot_data %>%
      filter(`Genus (gtdb)` %in% genus_levels_full)

    # y轴（基因）顺序：按全局总数降序
    gene_levels <- subdata %>%
      group_by(gene) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      filter(gene %in% top_genes) %>%
      arrange(desc(total)) %>%
      pull(gene)
    plot_data$gene <- factor(plot_data$gene, levels = gene_levels)
    
    # x轴Genus按出现总和排序（最多100）
    genus_levels <- plot_data %>%
      group_by(`Genus (gtdb)`) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      arrange(desc(total)) %>%
      pull(`Genus (gtdb)`)
    plot_data$`Genus (gtdb)` <- factor(plot_data$`Genus (gtdb)`, levels = genus_levels)
    
    # Host 为 facet
    n_genus <- length(unique(plot_data$`Genus (gtdb)`))
    n_gene <- length(gene_levels)
    n_host <- length(unique(plot_data$Host))
    fig_width <- 15
    fig_height <- 30

    # 安全文件名
    safe_contig <- gsub("[^A-Za-z0-9]", "_", contig_type)
    filename <- paste0(gene_name, "_", safe_contig, "_heatmap_by_host_TOP50.pdf")
    filepath <- file.path('Result/NCBI_4395_Batch/04_Gene_Taxonomy', filename)

    # 画图。色带从白色到#6566aa
    p <- ggplot(plot_data, aes(x = `Genus (gtdb)`, y = gene, fill = count)) +
      geom_tile(color = "white", size = 0.12) +
      scale_fill_gradient(
        low = "white", high = "#6566aa",
        name = "Count",
        trans = "log1p"
      ) +
      facet_wrap(~Host, nrow = n_host, scales = "fixed") +
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
        title = paste0(gene_name, " Top50 by Genus (Top100 Genus, Contig_Type: ", contig_type, ")"),
        x = "Genus (gtdb)",
        y = gene_name
      )
    
    ggsave(filepath, plot = p, width = fig_width, height = fig_height, limitsize = FALSE)
    
    # 保存当前Contig_Type的明细csv（只保留top50 gene和前100的Genus (gtdb)的行）
    csv_filename <- paste0(gene_name, "_", safe_contig, "_counts_by_host_TOP50.csv")
    csv_filepath <- file.path('Result/NCBI_4395_Batch/04_Gene_Taxonomy', csv_filename)
    write_csv(plot_data, csv_filepath)
  }
}

# 作图
create_heatmaps_by_contigtype(data, 'Defense_Subtype', 'Defense_Subtype')
create_heatmaps_by_contigtype(data, 'AntiDS_Type', 'AntiDS_Type')
create_heatmaps_by_contigtype(data, 'AMR_Type', 'AMR_Type')

cat("\nAll heatmaps completed!\n")

