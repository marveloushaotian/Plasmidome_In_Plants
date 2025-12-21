process CCTYPER_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/cctyper/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(gene_fna)

    output:
    tuple val(sample_id), path("${sample_id}_cas_operons.tab"), emit: cas_operons, optional: true
    tuple val(sample_id), path("${sample_id}_*"), emit: all_files
    path "${sample_id}_cctyper.log", emit: log

    script:
    """
    set -euo pipefail

    # Initialize log
    echo "[INFO] CCTyper analysis for ${sample_id}" > "${sample_id}_cctyper.log"
    echo "[INFO] Timestamp: \$(date)" >> "${sample_id}_cctyper.log"
    echo "[INFO] Input file: ${gene_fna}" >> "${sample_id}_cctyper.log"

    # Validate input
    if [ ! -f "${gene_fna}" ]; then
        echo "[ERROR] Input file does not exist" >> "${sample_id}_cctyper.log"
        touch "${sample_id}_cas_operons.tab"
        exit 0
    fi

    file_size=\$(stat -c%s "${gene_fna}" 2>/dev/null || echo 0)
    seq_count=\$(grep -c "^>" "${gene_fna}" 2>/dev/null || echo 0)

    echo "[INFO] File size: \$file_size bytes" >> "${sample_id}_cctyper.log"
    echo "[INFO] Sequences: \$seq_count" >> "${sample_id}_cctyper.log"

    # Handle empty input
    if [ \$file_size -eq 0 ] || [ \$seq_count -eq 0 ]; then
        echo "[WARN] Empty or invalid input file" >> "${sample_id}_cctyper.log"
        touch "${sample_id}_cas_operons.tab"
        exit 0
    fi

    # Check CCTyper availability
    if ! command -v cctyper &> /dev/null; then
        echo "[ERROR] CCTyper not found" >> "${sample_id}_cctyper.log"
        touch "${sample_id}_cas_operons.tab"
        exit 0
    fi

    # Run CCTyper directly in current directory
    echo "[INFO] Running CCTyper..." >> "${sample_id}_cctyper.log"

    # CCTyper expects output directory name, will create it
    # Use a unique temp name to avoid conflicts
    temp_output="cctyper_temp_\$(date +%s)_\$\$"

    # Run CCTyper with Prodigal meta mode and multi-threading
    timeout 1800 cctyper "${gene_fna}" "\$temp_output" --prodigal meta --threads ${task.cpus} > cctyper.stdout 2> cctyper.stderr || true

    # Log outputs
    if [ -f cctyper.stdout ]; then
        echo "[DEBUG] STDOUT:" >> "${sample_id}_cctyper.log"
        cat cctyper.stdout >> "${sample_id}_cctyper.log"
    fi
    if [ -f cctyper.stderr ]; then
        echo "[DEBUG] STDERR:" >> "${sample_id}_cctyper.log"
        cat cctyper.stderr >> "${sample_id}_cctyper.log"
    fi

    # Move all CCTyper output files to current directory WITH SAMPLE_ID PREFIX
    if [ -d "\$temp_output" ]; then
        echo "[INFO] CCTyper completed, collecting outputs..." >> "${sample_id}_cctyper.log"

        # Count and list output files
        file_count=\$(find "\$temp_output" -type f | wc -l)
        echo "[INFO] Found \$file_count output files" >> "${sample_id}_cctyper.log"

        # Move all files to current directory with sample_id prefix
        for file in "\$temp_output"/*; do
            if [ -f "\$file" ]; then
                filename=\$(basename "\$file")
                mv "\$file" "./${sample_id}_\$filename"
                echo "[INFO] Output file: ${sample_id}_\$filename" >> "${sample_id}_cctyper.log"
            fi
        done

        # Clean up temp directory
        rm -rf "\$temp_output"

        # Check if main result file exists
        if [ -f "${sample_id}_cas_operons.tab" ]; then
            result_count=\$(tail -n +2 "${sample_id}_cas_operons.tab" 2>/dev/null | wc -l || echo 0)
            echo "[INFO] Found \$result_count CRISPR-Cas systems in cas_operons.tab" >> "${sample_id}_cctyper.log"
        elif [ -f "${sample_id}_cas_operons_putative.tab" ]; then
            putative_count=\$(tail -n +2 "${sample_id}_cas_operons_putative.tab" 2>/dev/null | wc -l || echo 0)
            echo "[INFO] Found \$putative_count putative CRISPR-Cas systems in cas_operons_putative.tab" >> "${sample_id}_cctyper.log"
            echo "[INFO] No confirmed CRISPR-Cas systems found, only putative results available" >> "${sample_id}_cctyper.log"
            # Create empty cas_operons.tab file
            touch "${sample_id}_cas_operons.tab"
        else
            echo "[WARN] No CRISPR-Cas systems found" >> "${sample_id}_cctyper.log"
            touch "${sample_id}_cas_operons.tab"
        fi
    else
        echo "[WARN] CCTyper output directory not found, analysis may have failed" >> "${sample_id}_cctyper.log"
        touch "${sample_id}_cas_operons.tab"
    fi

    # Cleanup
    rm -f cctyper.stdout cctyper.stderr

    echo "[INFO] CCTyper analysis completed" >> "${sample_id}_cctyper.log"
    """
}
