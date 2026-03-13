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
    
    # Process all .tsv files (supports both CheckM1 and CheckM2 formats)
    for qa_file in ${qa_files}; do
      if [[ "\$qa_file" == *.tsv ]]; then
        # Detect file format: CheckM1 (Bin Id in header) or CheckM2 (Name in header)
        # Check first 2 lines to handle both merged files (with separator line) and single files
        first_line=\$(head -1 "\$qa_file")
        second_line=\$(head -2 "\$qa_file" | tail -1)
        
        if echo "\$first_line" | grep -q "Bin Id"; then
            # CheckM1 format: Completeness in column 7, Contamination in column 8
            comp_col=7
            cont_col=8
            echo "[INFO] Detected CheckM1 format in \$qa_file"
        elif echo "\$second_line" | grep -q "Bin Id"; then
            # CheckM1 merged format (first line is separator)
            comp_col=7
            cont_col=8
            echo "[INFO] Detected CheckM1 merged format in \$qa_file"
        elif echo "\$first_line" | grep -q "Name"; then
            # CheckM2 format: Completeness in column 2, Contamination in column 3
            comp_col=2
            cont_col=3
            echo "[INFO] Detected CheckM2 format in \$qa_file"
        else
            echo "[WARN] Unknown format in \$qa_file, skipping"
            continue
        fi
        
        awk -v comp_thresh=${params.completeness_threshold} -v cont_thresh=${params.contamination_threshold} \\
            -v comp_col=\$comp_col -v cont_col=\$cont_col '
          NR > 1 && NF > 2 && \$1 !~ /^[-]+\$/ && \$1 != "Bin" && \$1 != "Name" && \$1 !~ /^[[:space:]]*\$/ {
            genome=\$1
            # Use dynamic column access for awk
            for (i=1; i<=NF; i++) {
              if (i == comp_col) completeness = \$i
              if (i == cont_col) contamination = \$i
            }
            # Debug output
            printf "DEBUG: %s comp=%s cont=%s\\n", genome, completeness, contamination
            # Force numeric conversion
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
