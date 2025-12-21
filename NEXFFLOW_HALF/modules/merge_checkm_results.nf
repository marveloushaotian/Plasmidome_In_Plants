process MERGE_CHECKM_RESULTS {
    tag "merge-checkm"
    publishDir "${params.outdir}/checkm", mode: 'copy', overwrite: true
    
    input:
    path qa_all_files
    path qa_chro_files
    
    output:
    path "qa_all_merged.tsv", emit: qa_all
    path "qa_chro_merged.tsv", emit: qa_chro
    path "checkm_summary.txt", emit: summary
    
    script:
    """
    set -euo pipefail
    
    # 统计样本数量
    num_samples=\$(echo ${qa_all_files} | wc -w)
    echo "[INFO] Merging CheckM results from \$num_samples samples..."
    
    # 合并所有样本的 qa_all 结果
    # 第一个文件：保留表头
    first_all=\$(echo ${qa_all_files} | awk '{print \$1}')
    head -3 "\$first_all" > qa_all_merged.tsv
    
    # 所有文件：追加数据行（跳过前3行表头）
    for f in ${qa_all_files}; do
        tail -n +4 "\$f" >> qa_all_merged.tsv
    done
    
    # 合并所有样本的 qa_chro 结果
    first_chro=\$(echo ${qa_chro_files} | awk '{print \$1}')
    
    # 检查是否有有效的染色体结果
    if grep -q "Bin Id" "\$first_chro" 2>/dev/null; then
        head -3 "\$first_chro" > qa_chro_merged.tsv
        
        for f in ${qa_chro_files}; do
            # 跳过空文件或只有注释的文件
            if grep -q "Bin Id" "\$f" 2>/dev/null; then
                tail -n +4 "\$f" >> qa_chro_merged.tsv
            fi
        done
    else
        echo "# No valid chromosome sequences found in any sample" > qa_chro_merged.tsv
    fi
    
    # 生成汇总报告
    cat > checkm_summary.txt <<EOF
=== CheckM Results Summary ===
Date: \$(date)
Total samples processed: \$num_samples

Files merged:
- qa_all_merged.tsv: Quality assessment for all sequences
- qa_chro_merged.tsv: Quality assessment for chromosome-only sequences

Sample list:
EOF
    
    # 列出所有处理的样本
    for f in ${qa_all_files}; do
        basename "\$f" | sed 's/_qa_all.tsv//'
    done >> checkm_summary.txt
    
    echo "" >> checkm_summary.txt
    echo "Merge completed successfully." >> checkm_summary.txt
    """
}

