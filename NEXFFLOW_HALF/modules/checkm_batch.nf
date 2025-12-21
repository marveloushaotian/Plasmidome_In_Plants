process CHECKM_BATCH {
    tag "checkm-batch"
    publishDir "${params.outdir}/checkm", mode: 'copy', overwrite: true
    
    input:
    path bins_all_dir
    path bins_chro_dir
    
    output:
    path "qa_all.tsv",     emit: qa_all        // 改名避免冲突
    path "lineage_all.ms", emit: lineage_all   // 改名避免冲突
    path "qa_chro.tsv",    emit: qa_chro       // 改名避免冲突
    path "lineage_chro.ms",emit: lineage_chro  // 改名避免冲突
    
    script:
    """
    set -euo pipefail
    mkdir -p checkm_all checkm_chro
    
    checkm lineage_wf -x fasta "${bins_all_dir}" checkm_all -t ${task.cpus}
    checkm qa checkm_all/lineage.ms checkm_all -o 2 -f checkm_all/qa.tsv
    
    if compgen -G "${bins_chro_dir}/*.fasta" > /dev/null; then
      checkm lineage_wf -x fasta "${bins_chro_dir}" checkm_chro -t ${task.cpus}
      checkm qa checkm_chro/lineage.ms checkm_chro -o 2 -f checkm_chro/qa.tsv
    else
      echo "No chromosome-only FASTAs found" > checkm_chro/qa.tsv
    fi
    
    # 复制到不冲突的文件名
    cp checkm_all/qa.tsv qa_all.tsv
    cp checkm_all/lineage.ms lineage_all.ms
    cp checkm_chro/qa.tsv qa_chro.tsv
    cp checkm_chro/lineage.ms lineage_chro.ms
    """
}
