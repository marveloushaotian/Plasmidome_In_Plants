# Defense Systems Name Mapping File

## 文件位置
将Excel文件命名为 `Defense_Systems_Name_List.xlsx` 并放置在项目根目录的 `data/` 文件夹中。

## 文件格式

Excel文件应包含以下列：

| 列名 | 说明 | 示例 |
|------|------|------|
| deffinder_type | DefenseFinder 的类型名称 | CRISPR-Cas |
| deffinder_subtype | DefenseFinder 的子类型名称 | I-E |
| padloc_subtype | PADLOC 的子类型名称 | CRISPR-Cas_subtype_I-E |
| cctyper_subtype | CCTyper 的子类型名称 | I-E |
| unified type | 统一后的类型名称 | CRISPR-Cas |
| unified subtype | 统一后的子类型名称 | I-E |
| keep | 是否保留该系统 (TRUE/FALSE) | TRUE |

## 示例数据

```
deffinder_type,deffinder_subtype,padloc_subtype,cctyper_subtype,unified type,unified subtype,keep
CRISPR-Cas,I-E,CRISPR-Cas_subtype_I-E,I-E,CRISPR-Cas,I-E,TRUE
RM,I,RM_Type_I,NA,Restriction-Modification,Type-I,TRUE
Abi,AbiEii,AbiEii,NA,Abi,AbiEii,TRUE
DMS,DMS,DMS,NA,DMS,DMS,FALSE
```

## 说明

1. **统一命名**: 三个工具（DefenseFinder, PADLOC, CCTyper）对同一防御系统可能使用不同的名称。此映射文件用于统一命名。

2. **优先级**: 当多个工具检测到同一系统时，合并脚本会按照以下优先级选择：
   - CCTyper > DefenseFinder > PADLOC

3. **DMS和DRT**:
   - DMS (Defense system with Multiple Subsystems)
   - DRT (Defense-Related Terms)
   - 这些通常是不太特异的注释，在有其他选项时会被忽略

4. **keep 列**: 设置为 FALSE 的系统将从最终结果中过滤掉

## 创建方法

1. 收集三个工具的所有输出类型/子类型
2. 手动映射到统一的命名系统
3. 标记不需要的系统（keep=FALSE）
4. 保存为 Excel (.xlsx) 格式

## 注意事项

- 确保所有三个工具使用的名称都被映射
- NA 表示该工具不使用此分类
- 列名必须完全匹配（区分大小写）
