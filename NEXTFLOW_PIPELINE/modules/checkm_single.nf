process CHECKM_SINGLE {
    tag "${sample_id}"
    publishDir "${params.outdir}/checkm_single/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(labeled_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_qa_all.tsv"), emit: qa_all
    tuple val(sample_id), path("${sample_id}_qa_chro.tsv"), emit: qa_chro, optional: true
    path "${sample_id}_lineage_all.ms", emit: lineage_all
    path "${sample_id}_lineage_chro.ms", emit: lineage_chro, optional: true
    
    script:
    """
    set -euo pipefail
    
    # 为每个样本创建单独的输入目录
    mkdir -p input_all input_chro
    mkdir -p checkm_all checkm_chro
    
    # 复制文件到输入目录（CheckM 需要目录作为输入）
    cp "${labeled_fasta}" "input_all/${sample_id}.fasta"
    
    # 提取染色体序列
    awk '/^>/{keep=(\$0 ~ /\\|chromosome/)} keep' "${labeled_fasta}" > "input_chro/${sample_id}.fasta" || true
    
    # 运行 CheckM - 所有序列
    checkm lineage_wf -x fasta input_all checkm_all -t ${task.cpus}
    checkm qa checkm_all/lineage.ms checkm_all -o 2 -f checkm_all/qa.tsv
    
    # 运行 CheckM - 只有染色体（如果存在）
    if [ -s "input_chro/${sample_id}.fasta" ]; then
        checkm lineage_wf -x fasta input_chro checkm_chro -t ${task.cpus}
        checkm qa checkm_chro/lineage.ms checkm_chro -o 2 -f checkm_chro/qa.tsv
    else
        echo "# No chromosome sequences found for ${sample_id}" > checkm_chro/qa.tsv
        echo "# No chromosome sequences" > checkm_chro/lineage.ms
    fi
    
    # 复制输出文件并重命名（包含样本 ID 以避免文件名冲突）
    cp checkm_all/qa.tsv "${sample_id}_qa_all.tsv"
    cp checkm_chro/qa.tsv "${sample_id}_qa_chro.tsv"
    cp checkm_all/lineage.ms "${sample_id}_lineage_all.ms"
    cp checkm_chro/lineage.ms "${sample_id}_lineage_chro.ms"
    """
}

