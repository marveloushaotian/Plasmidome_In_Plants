process GENOMAD {
    tag "$sample_id"
    publishDir "${params.outdir}/genomad", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fasta)
    path genomad_db
    
    output:
    tuple val(sample_id), path("${sample_id}_genomad"), emit: summary
    path "${sample_id}_genomad/**/*.tsv", emit: tsv_files
    
    script:
    """
    genomad end-to-end \
        --cleanup \
        --splits 2 \
        --threads ${task.cpus} \
        ${fasta} \
        ${sample_id}_genomad \
        ${genomad_db}
    """
}
