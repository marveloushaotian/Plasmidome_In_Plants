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
    
    mkdir -p input_all input_chro
    mkdir -p checkm_all checkm_chro
    
    cp "${labeled_fasta}" "input_all/${sample_id}.fasta"
    
    awk '/^>/{keep=(\$0 ~ /\\|chromosome/)} keep' "${labeled_fasta}" > "input_chro/${sample_id}.fasta" || true
    
    checkm lineage_wf -x fasta input_all checkm_all -t ${task.cpus}
    checkm qa checkm_all/lineage.ms checkm_all -o 2 -f checkm_all/qa.tsv
    
    if [ -s "input_chro/${sample_id}.fasta" ]; then
        checkm lineage_wf -x fasta input_chro checkm_chro -t ${task.cpus}
        checkm qa checkm_chro/lineage.ms checkm_chro -o 2 -f checkm_chro/qa.tsv
    else
        echo "# No chromosome sequences found for ${sample_id}" > checkm_chro/qa.tsv
        echo "# No chromosome sequences" > checkm_chro/lineage.ms
    fi
    
    cp checkm_all/qa.tsv "${sample_id}_qa_all.tsv"
    cp checkm_chro/qa.tsv "${sample_id}_qa_chro.tsv"
    cp checkm_all/lineage.ms "${sample_id}_lineage_all.ms"
    cp checkm_chro/lineage.ms "${sample_id}_lineage_chro.ms"
    """
}

