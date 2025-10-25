process CCTYPER_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/cctyper/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(gene_fna)

    output:
    tuple val(sample_id), path("cas_operons.tab"), emit: cas_operons, optional: true
    tuple val(sample_id), path("*"), emit: all_files
    path "cctyper.log", emit: log

    script:
    """
    set -euo pipefail

    # Initialize log
    echo "[INFO] CCTyper analysis for ${sample_id}" > cctyper.log
    echo "[INFO] Timestamp: \$(date)" >> cctyper.log
    echo "[INFO] Input file: ${gene_fna}" >> cctyper.log

    # Validate input
    if [ ! -f "${gene_fna}" ]; then
        echo "[ERROR] Input file does not exist" >> cctyper.log
        touch cas_operons.tab
        exit 0
    fi

    file_size=\$(stat -c%s "${gene_fna}" 2>/dev/null || echo 0)
    seq_count=\$(grep -c "^>" "${gene_fna}" 2>/dev/null || echo 0)

    echo "[INFO] File size: \$file_size bytes" >> cctyper.log
    echo "[INFO] Sequences: \$seq_count" >> cctyper.log

    # Handle empty input
    if [ \$file_size -eq 0 ] || [ \$seq_count -eq 0 ]; then
        echo "[WARN] Empty or invalid input file" >> cctyper.log
        touch cas_operons.tab
        exit 0
    fi

    # Check CCTyper availability
    if ! command -v cctyper &> /dev/null; then
        echo "[ERROR] CCTyper not found" >> cctyper.log
        touch cas_operons.tab
        exit 0
    fi

    # Run CCTyper directly in current directory
    echo "[INFO] Running CCTyper..." >> cctyper.log

    # CCTyper expects output directory name, will create it
    # Use a unique temp name to avoid conflicts
    temp_output="cctyper_temp_\$(date +%s)_\$\$"

    # Run CCTyper
    timeout 1800 cctyper "${gene_fna}" "\$temp_output" > cctyper.stdout 2> cctyper.stderr || true

    # Log outputs
    if [ -f cctyper.stdout ]; then
        echo "[DEBUG] STDOUT:" >> cctyper.log
        cat cctyper.stdout >> cctyper.log
    fi
    if [ -f cctyper.stderr ]; then
        echo "[DEBUG] STDERR:" >> cctyper.log
        cat cctyper.stderr >> cctyper.log
    fi

    # Move all CCTyper output files to current directory
    if [ -d "\$temp_output" ]; then
        echo "[INFO] CCTyper completed, collecting outputs..." >> cctyper.log

        # Count and list output files
        file_count=\$(find "\$temp_output" -type f | wc -l)
        echo "[INFO] Found \$file_count output files" >> cctyper.log

        # Move all files to current directory
        for file in "\$temp_output"/*; do
            if [ -f "\$file" ]; then
                filename=\$(basename "\$file")
                mv "\$file" "./\$filename"
                echo "[INFO] Output file: \$filename" >> cctyper.log
            fi
        done

        # Clean up temp directory
        rm -rf "\$temp_output"

        # Check if main result file exists
        if [ -f "cas_operons.tab" ]; then
            result_count=\$(tail -n +2 "cas_operons.tab" 2>/dev/null | wc -l || echo 0)
            echo "[INFO] Found \$result_count CRISPR-Cas systems in cas_operons.tab" >> cctyper.log
        elif [ -f "cas_operons_putative.tab" ]; then
            echo "[INFO] Found cas_operons_putative.tab but not cas_operons.tab" >> cctyper.log
            # Use putative as main result if no confirmed operons
            cp "cas_operons_putative.tab" "cas_operons.tab"
        else
            echo "[WARN] No CRISPR-Cas systems found" >> cctyper.log
            touch cas_operons.tab
        fi
    else
        echo "[WARN] CCTyper output directory not found, analysis may have failed" >> cctyper.log
        touch cas_operons.tab
    fi

    # Cleanup
    rm -f cctyper.stdout cctyper.stderr

    echo "[INFO] CCTyper analysis completed" >> cctyper.log
    """
}
