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
    
    # Count samples
    num_samples=\$(echo ${qa_all_files} | wc -w)
    echo "[INFO] Merging CheckM results from \$num_samples samples..."
    
    # Merge all qa_all results
    # Detect format: CheckM1 (3-line header) or CheckM2 (1-line header)
    first_all=\$(echo ${qa_all_files} | awk '{print \$1}')
    
    if head -1 "\$first_all" | grep -q "^---"; then
        # CheckM1 format: 3-line header
        echo "[INFO] Detected CheckM1 format (3-line header)"
        head -3 "\$first_all" > qa_all_merged.tsv
        skip_lines=4
    elif head -1 "\$first_all" | grep -q "Name"; then
        # CheckM2 format: 1-line header
        echo "[INFO] Detected CheckM2 format (1-line header)"
        head -1 "\$first_all" > qa_all_merged.tsv
        skip_lines=2
    else
        echo "[WARN] Unknown format, using default (1-line header)"
        head -1 "\$first_all" > qa_all_merged.tsv
        skip_lines=2
    fi
    
    # Append data rows from all files (skip separator lines)
    for f in ${qa_all_files}; do
        tail -n +\$skip_lines "\$f" | grep -v "^---" | grep -v "^[[:space:]]*\$" >> qa_all_merged.tsv
    done
    
    # Merge all qa_chro results
    first_chro=\$(echo ${qa_chro_files} | awk '{print \$1}')
    
    # Check for valid chromosome results (CheckM1 or CheckM2 format)
    if head -1 "\$first_chro" | grep -q "^---"; then
        # CheckM1 format
        head -3 "\$first_chro" > qa_chro_merged.tsv
        skip_lines_chro=4
    elif head -1 "\$first_chro" | grep -q "Name"; then
        # CheckM2 format
        head -1 "\$first_chro" > qa_chro_merged.tsv
        skip_lines_chro=2
    elif grep -q "^#" "\$first_chro" 2>/dev/null; then
        # Empty file or comments
        echo "# No valid chromosome sequences found in any sample" > qa_chro_merged.tsv
        skip_lines_chro=0
    else
        echo "# No valid chromosome sequences found in any sample" > qa_chro_merged.tsv
        skip_lines_chro=0
    fi
    
    # Append data (skip separator lines and empty lines)
    if [ "\$skip_lines_chro" -gt 0 ]; then
        for f in ${qa_chro_files}; do
            # Skip empty files or files with only comments
            if ! grep -q "^#" "\$f" 2>/dev/null && [ -s "\$f" ]; then
                tail -n +\$skip_lines_chro "\$f" | grep -v "^---" | grep -v "^[[:space:]]*\$" >> qa_chro_merged.tsv
            fi
        done
    fi
    
    # Generate summary report
    cat > checkm_summary.txt <<EOF
=== CheckM Results Summary ===
Date: \$(date)
Total samples processed: \$num_samples

Files merged:
- qa_all_merged.tsv: Quality assessment for all sequences
- qa_chro_merged.tsv: Quality assessment for chromosome-only sequences

Sample list:
EOF
    
    # List all processed samples
    for f in ${qa_all_files}; do
        basename "\$f" | sed 's/_qa_all.tsv//'
    done >> checkm_summary.txt
    
    echo "" >> checkm_summary.txt
    echo "Merge completed successfully." >> checkm_summary.txt
    """
}

