# Defense Systems Merger Module 使用说明

## 功能概述

这个模块整合了三个防御系统注释工具的结果：
- **PADLOC** - 原核防御系统定位器
- **DefenseFinder** - 防御系统查找器
- **CCTyper** - CRISPR-Cas 分型工具

通过统一命名和智能合并，生成一致且全面的防御系统注释。

---

## 工作流程

```
HQ Genomes
    ↓
├→ PADLOC_ANNOTATE
├→ DEFENSEFINDER_ANNOTATE
└→ CCTYPER_ANNOTATE
    ↓
Collect & Merge
    ↓
├→ Collect GFF files (from Prokka)
├→ Collect PADLOC results
├→ Collect DefenseFinder results
└→ Collect CCTyper results
    ↓
Unify & Categorize
    ↓
├→ Defense Systems
├→ PDC Systems
└→ Antidefense Systems
    ↓
Final Output
```

---

## 准备工作

### 1. 创建映射文件

映射文件用于统一三个工具的命名系统。

**创建 Excel 文件：** `data/Defense_Systems_Name_List.xlsx`

**必需的列：**
- `deffinder_type` - DefenseFinder 类型
- `deffinder_subtype` - DefenseFinder 子类型
- `padloc_subtype` - PADLOC 子类型
- `cctyper_subtype` - CCTyper 子类型
- `unified type` - 统一的类型名称
- `unified subtype` - 统一的子类型名称
- `keep` - 是否保留 (TRUE/FALSE)

**示例内容：**
```
deffinder_type,deffinder_subtype,padloc_subtype,cctyper_subtype,unified type,unified subtype,keep
CRISPR-Cas,I-E,CRISPR-Cas_subtype_I-E,I-E,CRISPR-Cas,I-E,TRUE
RM,I,RM_Type_I,NA,Restriction-Modification,Type-I,TRUE
Abi,AbiEii,AbiEii,NA,Abi,AbiEii,TRUE
```

参考：`data/Defense_Systems_Mapping_Template.md`

### 2. 确保依赖项已运行

Defense merger 需要以下模块的输出：
- ✅ PROKKA_ANNOTATE (GFF 文件)
- ✅ PADLOC_ANNOTATE (防御系统)
- ✅ DEFENSEFINDER_ANNOTATE (防御系统)
- ✅ CCTYPER_ANNOTATE (CRISPR-Cas 系统)

这些会自动在 pipeline 中运行。

---

## 使用方法

### 基本用法

```bash
nextflow run main.nf \
  --input_mode labeled \
  --project_id "YOUR_PROJECT" \
  --run_defense_merger true \
  --defense_mapping_file "data/Defense_Systems_Name_List.xlsx"
```

### 完整示例

```bash
nextflow run main.nf \
  -profile standard \
  --input_mode labeled \
  --labeled_fastas "PREPARED_FROM_TK_ALL_LAEBLED" \
  --existing_qa "qa.tsv" \
  --project_id "Plasmidome_Defense_Analysis" \
  --run_defense_merger true \
  --defense_mapping_file "data/Defense_Systems_Name_List.xlsx" \
  --outdir "Results/Defense_Analysis"
```

---

## 输出文件

### 主输出目录
```
Results/YOUR_PROJECT/defense_merger/
```

### 文件说明

#### 1. **中间文件** (`intermediate/`)
- `all_transformed_gff.tsv` - 合并的 GFF 注释
- `final_padloc.tsv` - 合并的 PADLOC 结果
- `final_defensefinder.tsv` - 合并的 DefenseFinder 结果
- `final_cctyper.tsv` - 合并的 CCTyper 结果

#### 2. **合并结果**
- `all.merged` - 初步合并结果（包含所有工具的列）
- `all.merged.unified` - 统一命名后的完整结果

#### 3. **最终输出** ⭐
- `defense_results.tsv` - **主要结果文件**
  - 每个防御系统一行
  - 包含位置、类型、子类型、基因数量等信息

#### 4. **统计报告**
- `defense_summary_report.txt` - 文本格式摘要
- `defense_statistics.tsv` - 表格格式统计
- `merge_stats.txt` - 合并统计

---

## defense_results.tsv 列说明

| 列名 | 说明 |
|------|------|
| `Acc` | 基因组/contig ID |
| `Type` | 防御系统类型 |
| `Unified_subtype_name` | 统一的子类型名称 |
| `PDC_type` | PDC 系统类型 (如果有) |
| `PDC_unified_subtype_name` | PDC 子类型 |
| `Antidefense_type` | 反防御系统类型 (如果有) |
| `Antidefense_unified_subtype_name` | 反防御子类型 |
| `Start` | 起始位置 |
| `End` | 结束位置 |
| `sys_length` | 系统长度 (bp) |
| `gene_length` | 基因总长度 (bp) |
| `defense_num_genes` | 防御基因数量 |

---

## 合并逻辑

### 1. 命名统一
通过映射文件将三个工具的命名统一为标准名称。

### 2. 冲突解决
当多个工具检测到同一系统时，按优先级选择：
1. **CCTyper** (最高优先级 - CRISPR-Cas 专家)
2. **DefenseFinder** (中优先级 - 综合防御系统)
3. **PADLOC** (基础优先级 - 广泛检测)

### 3. 系统分类
- **Defense Systems** - 标准防御系统
- **PDC Systems** - 以 "PDC-" 开头的系统
- **Antidefense Systems** - DefenseFinder 标记为 "Antidefense" 的系统

### 4. 片段化处理
DefenseFinder 结果会被分割，如果蛋白质之间间隔 >5 个基因，则视为不同系统。

---

## 故障排除

### 错误：找不到映射文件
```
ERROR: Cannot find defense mapping file
```
**解决**：确保映射文件路径正确，使用绝对路径或相对于工作目录的路径。

### 错误：Python pandas 模块缺失
```
ModuleNotFoundError: No module named 'pandas'
```
**解决**：在 base conda 环境中安装依赖：
```bash
conda install -c conda-forge pandas openpyxl
```

### 警告：某些工具结果为空
```
[WARN] No PADLOC results found
```
**说明**：正常情况，说明该工具没有检测到相应的防御系统。

### 结果数量少于预期
**检查**：
1. 映射文件中的 `keep` 列是否设置为 TRUE
2. 是否有样本没有通过 HQ 筛选
3. 查看 `merge_stats.txt` 了解过滤统计

---

## 性能考虑

### 资源需求
- **CPU**: 4 cores (合并步骤)
- **内存**: 32 GB (合并步骤)
- **时间**: 通常 10-30 分钟 (取决于样本数量)

### 优化建议
1. 对于 >1000 个样本，考虑增加内存到 64 GB
2. 确保有足够的磁盘空间存储中间文件

---

## 示例输出

### defense_summary_report.txt
```
Defense Systems Analysis Summary
=================================

Total defense systems identified: 1523

By Type:
--------
   450 CRISPR-Cas
   312 Restriction-Modification
   198 Abi
   145 Toxin-Antitoxin
   ...

PDC Systems:
------------
Total PDC systems: 45

Antidefense Systems:
--------------------
Total antidefense systems: 23
```

---

## 参考资料

- **PADLOC**: https://github.com/padlocbio/padloc
- **DefenseFinder**: https://github.com/mdmparis/defense-finder
- **CCTyper**: https://github.com/Russel88/CCTyper

---

## 常见问题

**Q: 是否必须运行所有三个工具？**
A: 建议运行所有三个以获得最全面的结果，但模块可以处理部分缺失的情况。

**Q: 可以自定义优先级吗？**
A: 优先级硬编码在 Python 脚本中，修改需要编辑 `bin/merge_defense_systems.py`。

**Q: 映射文件可以用 CSV 格式吗？**
A: 目前只支持 Excel (.xlsx) 格式。

**Q: 如何只运行 merger 而不运行整个 pipeline？**
A: 暂不支持，merger 需要依赖其他模块的输出。

---

## 更新日志

- **v1.0** (2025-10-25): 初始版本，支持基本合并功能
