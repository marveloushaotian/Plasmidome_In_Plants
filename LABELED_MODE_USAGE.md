# 使用 Labeled 模式运行 Pipeline

## 概述

Pipeline 现在支持两种输入模式：

1. **RAW 模式**（原始模式）：从 FASTQ 原始测序数据开始，运行完整的 pipeline
2. **LABELED 模式**（新增）：从已经 labeled 的 FASTA 文件和 CheckM qa.tsv 开始

## Labeled 模式说明

当您已经有以下文件时，可以使用 labeled 模式：
- 已经通过 GeNomad 标记过的 FASTA 文件（位于某个目录下）
- CheckM 质量评估文件 (qa.tsv)

### 跳过的步骤

在 labeled 模式下，以下步骤会被跳过：
- Assembly（组装）
- Quality Filter（质量过滤）
- GeNomad（病毒/质粒检测）
- Label Contigs（contig 标记）
- CheckM（质量评估）

### 保留的步骤

所有下游分析保持不变：
- Filter HQ Genomes（高质量基因组筛选）
- GTDB-Tk（分类）
- MOB-suite（质粒分析）
- Prokka（注释）
- Prodigal（基因预测）
- PADLOC（防御系统）
- DefenseFinder（防御系统）
- CCTyper（CRISPR-Cas 系统）
- AMRFinder（耐药基因）
- EggNOG（功能注释）
- 系统发育树构建

### 质量筛选阈值

所有质量筛选阈值保持不变（在 main.nf 中定义）：
- Completeness (完整度): ≥ 90%
- Contamination (污染度): ≤ 5%

## 使用方法

### 1. 准备输入文件

确保您有以下文件结构：

```
Plasmidome_In_Plants/
├── PREPARED_FROM_TK_ALL_LAEBLED/    # 包含所有 labeled FASTA 文件的目录
│   ├── sample1.fasta
│   ├── sample2.fasta
│   ├── sample3.fasta
│   └── ...
├── qa.tsv                           # CheckM 质量评估文件
└── main.nf                          # Pipeline 主文件
```

**FASTA 文件命名要求：**
- 文件扩展名必须是 `.fasta`
- 文件名（去掉 .fasta 后缀）将作为样本 ID
- 例如：`sample1.fasta` → 样本 ID 为 `sample1`

**qa.tsv 格式要求：**
- 这应该是 CheckM 输出的标准 qa.tsv 文件
- 包含表头（前3行）
- 第1列：基因组/样本 ID
- 第7列：Completeness（完整度）
- 第8列：Contamination（污染度）

### 2. 运行 Pipeline

#### 基本用法（使用默认路径）

```bash
nextflow run main.nf \
  --input_mode labeled \
  --project_id "YOUR_PROJECT_NAME"
```

#### 自定义输入路径

```bash
nextflow run main.nf \
  --input_mode labeled \
  --labeled_fastas "path/to/your/labeled_fastas" \
  --existing_qa "path/to/your/qa.tsv" \
  --project_id "YOUR_PROJECT_NAME"
```

#### 完整示例（带所有重要参数）

```bash
nextflow run main.nf \
  --input_mode labeled \
  --labeled_fastas "PREPARED_FROM_TK_ALL_LAEBLED" \
  --existing_qa "qa.tsv" \
  --project_id "Plasmidome_Analysis" \
  --completeness_threshold 90.0 \
  --contamination_threshold 5.0 \
  --outdir "Results/Plasmidome_Analysis"
```

### 3. 运行系统发育树分析（可选）

#### 构建所有样本的整体系统发育树

```bash
nextflow run main.nf \
  --input_mode labeled \
  --project_id "YOUR_PROJECT_NAME" \
  --run_overall_tree true
```

#### 构建分组系统发育树

如果您想为不同组别的样本分别构建系统发育树：

```bash
nextflow run main.nf \
  --input_mode labeled \
  --project_id "YOUR_PROJECT_NAME" \
  --run_grouped_trees true \
  --sample_groups "sample_groups.tsv"
```

**sample_groups.tsv 格式：**
```
Sample_ID    Group_Name
sample1      GroupA
sample2      GroupA
sample3      GroupB
sample4      GroupB
```

### 4. 其他可选参数

```bash
# 启用 EggNOG 功能注释
--run_eggnog true

# 选择 EggNOG 使用的蛋白序列来源（prodigal 或 prokka）
--eggnog_protein_source prodigal

# 自定义数据库路径
--gtdbtk_db "/path/to/gtdbtk/database"
--eggnog_db "/path/to/eggnog/database"

# 自定义 IQ-TREE bootstrap 参数
--bootstrap_reps 1000
--alrt_reps 1000
```

## 输出结果

输出将保存在 `Results/YOUR_PROJECT_NAME/` 目录下，包含：

```
Results/YOUR_PROJECT_NAME/
├── hq_genomes/              # 高质量基因组列表
│   ├── hq_genomes.txt
│   ├── lq_genomes.txt
│   └── quality_summary.txt
├── gtdbtk/                  # GTDB-Tk 分类结果
├── mobsuite_recon/          # MOB-suite 重组分析
├── mobsuite_typer/          # MOB-suite 分型
├── prokka/                  # Prokka 注释
├── prodigal/                # Prodigal 基因预测
├── padloc/                  # PADLOC 防御系统
├── defensefinder/           # DefenseFinder 结果
├── cctyper/                 # CCTyper CRISPR-Cas 分析
├── amrfinder/               # AMRFinder 耐药基因
├── eggnog/                  # EggNOG 功能注释（如果启用）
├── iqtree/                  # 系统发育树（如果启用）
└── pipeline_info/           # Pipeline 运行信息
```

## 切换回 RAW 模式

如果需要从原始 FASTQ 数据运行完整 pipeline：

```bash
nextflow run main.nf \
  --input_mode raw \
  --raw_reads "Raw_Reads" \
  --project_id "YOUR_PROJECT_NAME"
```

## 注意事项

1. **样本 ID 匹配**：确保 qa.tsv 中的样本 ID 与 FASTA 文件名（去掉 .fasta）一致
2. **文件格式**：FASTA 文件应该已经包含 GeNomad 的标记（如 `|chromosome`, `|plasmid` 等）
3. **内存和 CPU**：根据您的样本数量和服务器配置，可能需要在 `nextflow.config` 中调整资源分配
4. **数据库**：确保已经正确设置了 GTDB-Tk 和 EggNOG 数据库路径

## 故障排除

### 错误：找不到输入文件

确保：
- `PREPARED_FROM_TK_ALL_LAEBLED` 目录存在
- 目录中包含 `.fasta` 文件
- `qa.tsv` 文件存在且格式正确

### 错误：没有 HQ 基因组

检查：
- `qa.tsv` 文件格式是否正确
- 是否有样本满足质量阈值（Completeness ≥ 90%, Contamination ≤ 5%）
- 可以尝试降低阈值：`--completeness_threshold 50 --contamination_threshold 10`

### 错误：样本 ID 不匹配

确保 FASTA 文件名与 qa.tsv 中的样本 ID 一致。
