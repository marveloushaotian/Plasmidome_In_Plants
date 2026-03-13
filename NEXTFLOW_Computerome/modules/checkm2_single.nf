process CHECKM2_SINGLE {
    tag "${sample_id}"
    publishDir "${params.outdir}/checkm2_single/${sample_id}", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(labeled_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_qa_all.tsv"), emit: qa_all
    tuple val(sample_id), path("${sample_id}_qa_chro.tsv"), emit: qa_chro, optional: true
    
    script:
    """
    set -euo pipefail
    
    # Set CheckM2 database path
    export CHECKM2DB="${params.checkm2_db}"
    
    # Create separate input directories for each sample
    mkdir -p input_all input_chro
    
    # Copy file to input directory
    cp "${labeled_fasta}" "input_all/${sample_id}.fasta"
    
    # Extract chromosome sequences
    awk '/^>/{keep=(\$0 ~ /\\|chromosome/)} keep' "${labeled_fasta}" > "input_chro/${sample_id}.fasta" || true
    
    # Run CheckM2 - all sequences
    checkm2 predict \\
        --threads ${task.cpus} \\
        --input input_all \\
        --output-directory checkm2_all \\
        -x fasta
    
    # Copy CheckM2 raw output directly without any conversion
    cp checkm2_all/quality_report.tsv ${sample_id}_qa_all.tsv
    
    echo "[INFO] Copied CheckM2 output for ${sample_id} (all sequences)"
    
    # Run CheckM2 - chromosome only (if exists)
    if [ -s "input_chro/${sample_id}.fasta" ]; then
        checkm2 predict \\
            --threads ${task.cpus} \\
            --input input_chro \\
            --output-directory checkm2_chro \\
            -x fasta
        
        # Copy chromosome results directly
        if [ -f "checkm2_chro/quality_report.tsv" ]; then
            cp checkm2_chro/quality_report.tsv ${sample_id}_qa_chro.tsv
            echo "[INFO] Copied CheckM2 output for ${sample_id} (chromosome only)"
        else
            echo "# No chromosome sequences found for ${sample_id}" > ${sample_id}_qa_chro.tsv
            echo "[WARN] No valid chromosome data for ${sample_id}"
        fi
    else
        echo "# No chromosome sequences found for ${sample_id}" > ${sample_id}_qa_chro.tsv
    fi
    """
}

