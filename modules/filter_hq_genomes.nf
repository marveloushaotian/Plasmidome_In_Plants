process FILTER_HQ_GENOMES {
    tag "filter-hq"
    publishDir "${params.outdir}/hq_genomes", mode: 'copy', overwrite: true
    
    input:
    path qa_files
    
    output:
    path "hq_genomes.txt",      emit: hq_list
    path "lq_genomes.txt",      emit: lq_list
    path "quality_summary.txt", emit: summary
    
    script:
    """
    set -euo pipefail
    > hq_genomes.txt
    > lq_genomes.txt
    > quality_summary.txt
    
    # 只处理 qa_all 文件（支持 qa_all.tsv 和 qa_all_merged.tsv），避免重复
    for qa_file in ${qa_files}; do
      if [[ "\$qa_file" == *"qa_all"*".tsv" ]]; then
        awk -v comp_thresh=${params.completeness_threshold} -v cont_thresh=${params.contamination_threshold} '
          NF > 20 && \$1 !~ /^[-]+\$/ && \$1 != "Bin" && \$1 !~ /^[[:space:]]*\$/ {
            genome=\$1
            completeness=\$7
            contamination=\$8
            # 调试输出
            printf "DEBUG: %s comp=%s cont=%s\\n", genome, completeness, contamination
            # 强制数值转换
            comp_val = completeness + 0.0
            cont_val = contamination + 0.0
            printf "DEBUG: %s comp_val=%.2f cont_val=%.2f\\n", genome, comp_val, cont_val
            if (comp_val >= comp_thresh && cont_val <= cont_thresh) {
              print genome >> "hq_genomes.txt"
              printf "DEBUG: %s -> HQ\\n", genome
            } else {
              print genome >> "lq_genomes.txt"
              printf "DEBUG: %s -> LQ\\n", genome
            }
          }
        ' "\$qa_file"
      fi
    done
    
    echo "High quality genomes: \$(wc -l < hq_genomes.txt)" >> quality_summary.txt
    echo "Low  quality genomes: \$(wc -l < lq_genomes.txt)" >> quality_summary.txt
    """
}
