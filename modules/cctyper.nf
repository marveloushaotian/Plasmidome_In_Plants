process CCTYPER_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/cctyper/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(gene_fna)
    
    output:
    tuple val(sample_id), path("${sample_id}_cctyper_results.tsv"), emit: results
    tuple val(sample_id), path("${sample_id}_*"), emit: all_files, optional: true
    path "${sample_id}_cctyper.log", emit: log
    
    script:
    """
    set -euo pipefail
    
    # Initialize log with basic info
    echo "[INFO] CCTyper analysis for ${sample_id}" > "${sample_id}_cctyper.log"
    echo "[INFO] Timestamp: \$(date)" >> "${sample_id}_cctyper.log"
    echo "[INFO] Input file: ${gene_fna}" >> "${sample_id}_cctyper.log"
    
    # Quick input validation
    if [ ! -f "${gene_fna}" ]; then
        echo "[ERROR] Input file does not exist: ${gene_fna}" >> "${sample_id}_cctyper.log"
        echo "# CCTyper failed - input file missing" > "${sample_id}_cctyper_results.tsv"
        echo "Contig\\tOperon\\tOperon_Pos\\tCompleteness\\tComplete_Interference\\tCas_Genes\\tOther_Genes\\tCRISPR_Genes\\tCas_Subtype\\tCRISPR_Subtype\\tPrediction" >> "${sample_id}_cctyper_results.tsv"
        exit 0
    fi
    
    file_size=\$(stat -c%s "${gene_fna}" 2>/dev/null || echo 0)
    seq_count=\$(grep -c "^>" "${gene_fna}" 2>/dev/null || echo 0)
    
    echo "[INFO] File size: \$file_size bytes" >> "${sample_id}_cctyper.log"
    echo "[INFO] Sequences: \$seq_count" >> "${sample_id}_cctyper.log"
    
    # Handle empty or invalid input
    if [ \$file_size -eq 0 ] || [ \$seq_count -eq 0 ]; then
        echo "[WARN] Empty or invalid input file" >> "${sample_id}_cctyper.log"
        echo "# No genes available for ${sample_id}" > "${sample_id}_cctyper_results.tsv"
        echo "Contig\\tOperon\\tOperon_Pos\\tCompleteness\\tComplete_Interference\\tCas_Genes\\tOther_Genes\\tCRISPR_Genes\\tCas_Subtype\\tCRISPR_Subtype\\tPrediction" >> "${sample_id}_cctyper_results.tsv"
        exit 0
    fi
    
    # Validate input file size (seems suspicious if very small but has many sequences)
    if [ \$file_size -lt 1000 ] && [ \$seq_count -gt 1000 ]; then
        echo "[WARN] Suspicious input file: small size (\$file_size bytes) but many sequences (\$seq_count)" >> "${sample_id}_cctyper.log"
    fi
    
    # Check CCTyper availability
    if ! command -v cctyper &> /dev/null; then
        echo "[ERROR] CCTyper not found" >> "${sample_id}_cctyper.log"
        echo "# CCTyper not available" > "${sample_id}_cctyper_results.tsv"
        echo "Contig\\tOperon\\tOperon_Pos\\tCompleteness\\tComplete_Interference\\tCas_Genes\\tOther_Genes\\tCRISPR_Genes\\tCas_Subtype\\tCRISPR_Subtype\\tPrediction" >> "${sample_id}_cctyper_results.tsv"
        exit 0
    fi
    
    # Create unique output directory name using mktemp in work directory
    # Use a name that includes timestamp and process ID for uniqueness
    work_dir="cctyper_${sample_id}_\$(date +%s)_\$\$"
    
    echo "[INFO] Creating work directory: \$work_dir" >> "${sample_id}_cctyper.log"
    mkdir -p "\$work_dir"
    output_dir="\$work_dir/cctyper_output"
    
    echo "[INFO] Output directory: \$output_dir" >> "${sample_id}_cctyper.log"
    echo "[INFO] Directory created at: \$(date)" >> "${sample_id}_cctyper.log"
    
    # Verify directory doesn't exist (shouldn't, we just created parent)
    if [ -d "\$output_dir" ]; then
        echo "[WARN] Output directory already exists, removing..." >> "${sample_id}_cctyper.log"
        rm -rf "\$output_dir"
    fi
    
    file_count=\$(find "\$work_dir" -type f | wc -l)
    echo "[INFO] Files in work directory before cctyper: \$file_count" >> "${sample_id}_cctyper.log"
    
    # Run CCTyper with timeout - redirect STDERR to filter out "Directory already exists" message
    echo "[INFO] Running CCTyper..." >> "${sample_id}_cctyper.log"
    timeout 1800 cctyper "${gene_fna}" "\$output_dir" > cctyper.out 2> cctyper.err
    exit_code=\$?
    
    echo "[INFO] CCTyper exit code: \$exit_code" >> "${sample_id}_cctyper.log"
    
    # Check if cctyper actually failed or just complained about directory
    if grep -q "ERROR.*Directory.*already exists" cctyper.err 2>/dev/null; then
        echo "[WARN] CCTyper reported directory exists, but continuing..." >> "${sample_id}_cctyper.log"
        # If exit code is 0 and directory has files, consider it successful
        if [ \$exit_code -eq 0 ]; then
            echo "[INFO] CCTyper exit code is 0, treating as success despite directory warning" >> "${sample_id}_cctyper.log"
        fi
    fi
    
    # Log outputs
    if [ -f cctyper.out ]; then
        echo "[DEBUG] STDOUT:" >> "${sample_id}_cctyper.log"
        cat cctyper.out >> "${sample_id}_cctyper.log"
    fi
    if [ -f cctyper.err ]; then
        echo "[DEBUG] STDERR:" >> "${sample_id}_cctyper.log"
        cat cctyper.err >> "${sample_id}_cctyper.log"
    fi
    
    # Process results
    if [ \$exit_code -eq 0 ] && [ -d "\$output_dir" ]; then
        echo "[INFO] CCTyper completed successfully" >> "${sample_id}_cctyper.log"
        
        # Copy output files
        for file in "\$output_dir"/*; do
            if [ -f "\$file" ]; then
                filename=\$(basename "\$file")
                cp "\$file" "${sample_id}_\$filename"
                echo "[INFO] Output: \$filename" >> "${sample_id}_cctyper.log"
            fi
        done
        
        # Create main results file
        if [ -f "\$output_dir/cas_operons_putative.tab" ]; then
            cp "\$output_dir/cas_operons_putative.tab" "${sample_id}_cctyper_results.tsv"
        elif [ -f "\$output_dir/cas_operons.tab" ]; then
            cp "\$output_dir/cas_operons.tab" "${sample_id}_cctyper_results.tsv"
        else
            # CCTyper succeeded but found no CRISPR systems
            echo "# No CRISPR-Cas systems found by CCTyper" > "${sample_id}_cctyper_results.tsv"
            echo "Contig\\tOperon\\tOperon_Pos\\tCompleteness\\tComplete_Interference\\tCas_Genes\\tOther_Genes\\tCRISPR_Genes\\tCas_Subtype\\tCRISPR_Subtype\\tPrediction" >> "${sample_id}_cctyper_results.tsv"
        fi
    else
        echo "[WARN] CCTyper failed or timed out" >> "${sample_id}_cctyper.log"
        echo "# CCTyper analysis failed" > "${sample_id}_cctyper_results.tsv"
        echo "Contig\\tOperon\\tOperon_Pos\\tCompleteness\\tComplete_Interference\\tCas_Genes\\tOther_Genes\\tCRISPR_Genes\\tCas_Subtype\\tCRISPR_Subtype\\tPrediction" >> "${sample_id}_cctyper_results.tsv"
    fi
    
    # Report final stats
    result_lines=\$(grep -v "^#" "${sample_id}_cctyper_results.tsv" | tail -n +2 | wc -l 2>/dev/null || echo 0)
    echo "[INFO] Final result: \$result_lines CRISPR-Cas systems found" >> "${sample_id}_cctyper.log"
    
    # Cleanup
    rm -rf "\$work_dir" cctyper.out cctyper.err
    """
}

process CCTYPER_SUMMARY {
    tag "cctyper-summary"
    publishDir "${params.outdir}/cctyper_summary", mode: 'copy', overwrite: true
    
    input:
    path result_files
    
    output:
    path "cctyper_summary.tsv", emit: summary
    path "crispr_cas_statistics.txt", emit: statistics
    
    when:
    result_files.size() > 0
    
    script:
    """
    set -euo pipefail
    
    # Create summary
    echo "Sample_ID\\tTotal_CRISPR\\tType_I\\tType_II\\tType_III\\tType_IV\\tType_V\\tType_VI\\tOther\\tStatus" > cctyper_summary.tsv
    
    # Statistics file
    echo "CCTyper CRISPR-Cas Analysis Summary" > crispr_cas_statistics.txt
    echo "===================================" >> crispr_cas_statistics.txt
    echo "Generated: \$(date)" >> crispr_cas_statistics.txt
    echo "" >> crispr_cas_statistics.txt
    
    total_samples=0
    successful_samples=0
    failed_samples=0
    samples_with_crispr=0
    total_systems=0
    
    for result_file in ${result_files}; do
        if [ -f "\$result_file" ]; then
            sample_id=\$(basename "\$result_file" _cctyper_results.tsv)
            total_samples=\$((total_samples + 1))
            
            # Determine status
            if grep -q "^# CCTyper failed\\|^# CCTyper not available\\|^# No genes" "\$result_file"; then
                status="Failed"
                failed_samples=\$((failed_samples + 1))
                crispr_count=0
            else
                status="Success"
                successful_samples=\$((successful_samples + 1))
                crispr_count=\$(grep -v "^#" "\$result_file" | tail -n +2 | wc -l 2>/dev/null || echo 0)
            fi
            
            if [ \$crispr_count -gt 0 ]; then
                samples_with_crispr=\$((samples_with_crispr + 1))
                total_systems=\$((total_systems + crispr_count))
            fi
            
            # Simple type counting (can be improved based on actual output format)
            type_counts="0\\t0\\t0\\t0\\t0\\t0\\t\$crispr_count"
            
            echo "\$sample_id\\t\$crispr_count\\t\$type_counts\\t\$status" >> cctyper_summary.tsv
            
            echo "Sample \$sample_id: \$status (\$crispr_count systems)" >> crispr_cas_statistics.txt
        fi
    done
    
    echo "" >> crispr_cas_statistics.txt
    echo "Overall Statistics:" >> crispr_cas_statistics.txt
    echo "- Total samples: \$total_samples" >> crispr_cas_statistics.txt
    echo "- Successful: \$successful_samples" >> crispr_cas_statistics.txt
    echo "- Failed: \$failed_samples" >> crispr_cas_statistics.txt
    echo "- With CRISPR systems: \$samples_with_crispr" >> crispr_cas_statistics.txt
    echo "- Total CRISPR systems: \$total_systems" >> crispr_cas_statistics.txt
    
    if [ \$successful_samples -gt 0 ]; then
        avg=\$(echo "\$total_systems \$successful_samples" | awk '{printf "%.2f", \$1 / \$2}')
        echo "- Average per successful sample: \$avg" >> crispr_cas_statistics.txt
    fi
    """
}