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
    
    num_samples=\$(echo ${qa_all_files} | wc -w)
    echo "[INFO] Merging CheckM results from \$num_samples samples..."
    
    first_all=\$(echo ${qa_all_files} | awk '{print \$1}')
    head -3 "\$first_all" > qa_all_merged.tsv
    
    for f in ${qa_all_files}; do
        tail -n +4 "\$f" >> qa_all_merged.tsv
    done
    
    first_chro=\$(echo ${qa_chro_files} | awk '{print \$1}')
    
    if grep -q "Bin Id" "\$first_chro" 2>/dev/null; then
        head -3 "\$first_chro" > qa_chro_merged.tsv
        
        for f in ${qa_chro_files}; do
            if grep -q "Bin Id" "\$f" 2>/dev/null; then
                tail -n +4 "\$f" >> qa_chro_merged.tsv
            fi
        done
    else
        echo "# No valid chromosome sequences found in any sample" > qa_chro_merged.tsv
    fi
    
    cat > checkm_summary.txt <<EOF
=== CheckM Results Summary ===
Date: \$(date)
Total samples processed: \$num_samples

Files merged:
- qa_all_merged.tsv: Quality assessment for all sequences
- qa_chro_merged.tsv: Quality assessment for chromosome-only sequences

Sample list:
EOF
    
    for f in ${qa_all_files}; do
        basename "\$f" | sed 's/_qa_all.tsv//'
    done >> checkm_summary.txt
    
    echo "" >> checkm_summary.txt
    echo "Merge completed successfully." >> checkm_summary.txt
    """
}

