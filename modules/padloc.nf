process PADLOC_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/padloc/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(nucleotide_fna)
    
    output:
    tuple val(sample_id), path("${sample_id}_padloc.csv"), emit: results
    tuple val(sample_id), path("${sample_id}_*"), emit: all_files
    path "${sample_id}_padloc.log", emit: log
    
    script:
    """
    set -euo pipefail
    
    echo "[INFO] Starting PADLOC defense system annotation for ${sample_id}" > "${sample_id}_padloc.log"
    echo "[INFO] Timestamp: \$(date)" >> "${sample_id}_padloc.log"
    echo "[INFO] Working directory: \$(pwd)" >> "${sample_id}_padloc.log"
    echo "[INFO] Using nucleotide sequence mode" >> "${sample_id}_padloc.log"
    
    if [ ! -s "${nucleotide_fna}" ]; then
        echo "[WARN] Input nucleotide file ${nucleotide_fna} is empty or doesn't exist" >&2
        echo "[WARN] Input nucleotide file ${nucleotide_fna} is empty or doesn't exist" >> "${sample_id}_padloc.log"
        touch "${sample_id}_padloc.domtblout"
        exit 0
    fi
    
    if ! command -v padloc &> /dev/null; then
        echo "[ERROR] padloc not found in PATH" >&2
        echo "[ERROR] padloc not found in PATH" >> "${sample_id}_padloc.log"
        touch "${sample_id}_padloc.domtblout"
        exit 0
    fi
    
    echo "[INFO] Nucleotide file: ${nucleotide_fna}" >> "${sample_id}_padloc.log"
    echo "[INFO] Sample ID: ${sample_id}" >> "${sample_id}_padloc.log"
    
    echo "[INFO] PADLOC version:" >> "${sample_id}_padloc.log"
    padloc --version >> "${sample_id}_padloc.log" 2>&1 || echo "[WARN] Version check failed" >> "${sample_id}_padloc.log"
    
    echo "[INFO] Input file information:" >> "${sample_id}_padloc.log"
    echo "  File size: \$(stat -c%s '${nucleotide_fna}' 2>/dev/null || echo 'unknown') bytes" >> "${sample_id}_padloc.log"
    echo "  Sequences count: \$(grep -c '^>' '${nucleotide_fna}' 2>/dev/null || echo 'unknown')" >> "${sample_id}_padloc.log"
    
    mkdir -p padloc_output
    echo "[INFO] Created output directory: padloc_output" >> "${sample_id}_padloc.log"
    
    echo "[INFO] Running PADLOC command:" >> "${sample_id}_padloc.log"
    echo "padloc --fna ${nucleotide_fna} --outdir padloc_output --cpu ${task.cpus}" >> "${sample_id}_padloc.log"
    echo "[INFO] Starting PADLOC execution at: \$(date)" >> "${sample_id}_padloc.log"
    
    padloc --fna "${nucleotide_fna}" --outdir padloc_output --cpu ${task.cpus} > padloc_stdout.log 2> padloc_stderr.log
    exit_code=\$?
    
    echo "[INFO] PADLOC finished at: \$(date)" >> "${sample_id}_padloc.log"
    echo "[INFO] PADLOC exit code: \$exit_code" >> "${sample_id}_padloc.log"
    
    echo "[INFO] PADLOC stdout:" >> "${sample_id}_padloc.log"
    cat padloc_stdout.log >> "${sample_id}_padloc.log" 2>/dev/null || echo "No stdout captured" >> "${sample_id}_padloc.log"
    echo "[INFO] PADLOC stderr:" >> "${sample_id}_padloc.log"
    cat padloc_stderr.log >> "${sample_id}_padloc.log" 2>/dev/null || echo "No stderr captured" >> "${sample_id}_padloc.log"
    
    echo "[INFO] Checking for output files:" >> "${sample_id}_padloc.log"
    ls -la padloc_output/ >> "${sample_id}_padloc.log" 2>&1 || echo "No padloc_output directory found" >> "${sample_id}_padloc.log"
    
    csv_found=false
    csv_source=""
    
    if [ -f "padloc_output/${sample_id}_padloc.csv" ] && [ -s "padloc_output/${sample_id}_padloc.csv" ]; then
        csv_source="padloc_output/${sample_id}_padloc.csv"
        csv_found=true
    else
        csv_file=\$(find padloc_output/ -name "*_padloc.csv" -type f -size +0c | head -1)
        if [ -n "\$csv_file" ] && [ -f "\$csv_file" ]; then
            csv_source="\$csv_file"
            csv_found=true
        fi
    fi
    
    if [ "\$csv_found" = true ] && [ \$exit_code -eq 0 ]; then
        cp "\$csv_source" "${sample_id}_padloc.csv"
        
        output_size=\$(stat -c%s "${sample_id}_padloc.csv" 2>/dev/null || echo 0)
        systems_count=\$(tail -n +2 "${sample_id}_padloc.csv" | wc -l 2>/dev/null || echo 0)
        
        echo "[INFO] PADLOC completed successfully for ${sample_id}" >&2
        echo "[INFO] Output file size: \$output_size bytes" >> "${sample_id}_padloc.log"
        echo "[INFO] Defense systems found: \$systems_count" >> "${sample_id}_padloc.log"
        
        for file in padloc_output/*; do
            if [ -f "\$file" ]; then
                filename=\$(basename "\$file")
                if [[ "\$filename" != *"${sample_id}"* ]]; then
                    cp "\$file" "${sample_id}_\$filename"
                else
                    cp "\$file" "\$filename"
                fi
            fi
        done
        
    else
        echo "[WARN] PADLOC failed or no output found for ${sample_id}" >&2
        echo "[WARN] PADLOC exit code: \$exit_code" >> "${sample_id}_padloc.log"
        touch "${sample_id}_padloc.csv"
    fi
    
    rm -f padloc_stdout.log padloc_stderr.log
    
    echo "[INFO] Process completed at: \$(date)" >> "${sample_id}_padloc.log"
    """
}

process PADLOC_SUMMARY {
    tag "padloc-summary"
    publishDir "${params.outdir}/padloc_summary", mode: 'copy', overwrite: true
    
    input:
    path result_files
    
    output:
    path "padloc_summary.tsv", emit: summary
    path "defense_systems_statistics.txt", emit: statistics
    
    script:
    """
    set -euo pipefail
    
    echo "Sample_ID\\tTotal_Systems\\tCRISPR\\tRestriction_Modification\\tAbi\\tToxin_Antitoxin\\tOther\\tFile_Size_Bytes" > padloc_summary.tsv
    
    echo "PADLOC Defense Systems Summary (FNA Mode)" > defense_systems_statistics.txt
    echo "=========================================" >> defense_systems_statistics.txt
    echo "Generated on: \$(date)" >> defense_systems_statistics.txt
    echo "" >> defense_systems_statistics.txt
    
    total_samples=0
    samples_with_systems=0
    total_systems=0
    total_crispr=0
    total_rm=0
    total_abi=0
    total_ta=0
    total_other=0
    
    for result_file in ${result_files}; do
        if [ -f "\$result_file" ]; then
            sample_id=\$(basename "\$result_file" _padloc.csv)
            file_size=\$(stat -c%s "\$result_file" 2>/dev/null || echo 0)
            
            total_systems_count=0
            if [ -s "\$result_file" ]; then
                total_systems_count=\$(tail -n +2 "\$result_file" | wc -l 2>/dev/null || echo 0)
            fi
            
            crispr_count=0
            rm_count=0
            abi_count=0
            ta_count=0
            other_count=0
            
            if [ \$total_systems_count -gt 0 ]; then
                samples_with_systems=\$((samples_with_systems + 1))
                
                crispr_count=\$(tail -n +2 "\$result_file" | cut -d',' -f3 | grep -i "crispr\\|cas\\|type.*i\\|type.*ii\\|type.*iii\\|type.*iv\\|type.*v\\|type.*vi" | wc -l 2>/dev/null || echo 0)
                
                rm_count=\$(tail -n +2 "\$result_file" | cut -d',' -f3 | grep -i "restriction\\|modification\\|rm\\|hsdm\\|hsdr\\|hsds" | wc -l 2>/dev/null || echo 0)
                
                abi_count=\$(tail -n +2 "\$result_file" | cut -d',' -f3 | grep -i "abortive\\|abi\\|abort" | wc -l 2>/dev/null || echo 0)
                
                ta_count=\$(tail -n +2 "\$result_file" | cut -d',' -f3 | grep -i "toxin\\|antitoxin\\|mazf\\|rele\\|hipa\\|vapb\\|vapc" | wc -l 2>/dev/null || echo 0)
                
                classified_count=\$((crispr_count + rm_count + abi_count + ta_count))
                if [ \$classified_count -le \$total_systems_count ]; then
                    other_count=\$((total_systems_count - classified_count))
                else
                    other_count=0
                fi
            fi
            
            echo "\$sample_id\\t\$total_systems_count\\t\$crispr_count\\t\$rm_count\\t\$abi_count\\t\$ta_count\\t\$other_count\\t\$file_size" >> padloc_summary.tsv
            
            echo "Sample: \$sample_id" >> defense_systems_statistics.txt
            echo "  Total defense systems: \$total_systems_count" >> defense_systems_statistics.txt
            if [ \$total_systems_count -gt 0 ]; then
                echo "  System breakdown:" >> defense_systems_statistics.txt
                tail -n +2 "\$result_file" | cut -d',' -f3 | sort | uniq -c | sed 's/^/    /' >> defense_systems_statistics.txt
            else
                echo "  No defense systems found" >> defense_systems_statistics.txt
            fi
            echo "" >> defense_systems_statistics.txt
            
            total_samples=\$((total_samples + 1))
            total_systems=\$((total_systems + total_systems_count))
            total_crispr=\$((total_crispr + crispr_count))
            total_rm=\$((total_rm + rm_count))
            total_abi=\$((total_abi + abi_count))
            total_ta=\$((total_ta + ta_count))
            total_other=\$((total_other + other_count))
        fi
    done
    
    if [ \$total_samples -gt 0 ]; then
        avg_systems=\$(echo "scale=2; \$total_systems / \$total_samples" | bc -l 2>/dev/null || echo "0")
    else
        avg_systems=0
    fi
    
    echo "Overall Statistics:" >> defense_systems_statistics.txt
    echo "  Total samples: \$total_samples" >> defense_systems_statistics.txt
    echo "  Samples with defense systems: \$samples_with_systems" >> defense_systems_statistics.txt
    echo "  Total defense systems: \$total_systems" >> defense_systems_statistics.txt
    echo "  Average systems per sample: \$avg_systems" >> defense_systems_statistics.txt
    echo "  CRISPR systems: \$total_crispr" >> defense_systems_statistics.txt
    echo "  RM systems: \$total_rm" >> defense_systems_statistics.txt
    echo "  Abi systems: \$total_abi" >> defense_systems_statistics.txt
    echo "  TA systems: \$total_ta" >> defense_systems_statistics.txt
    echo "  Other systems: \$total_other" >> defense_systems_statistics.txt
    echo "" >> defense_systems_statistics.txt
    echo "Note: FNA mode analysis provides system-level identification." >> defense_systems_statistics.txt
    echo "Individual systems may contain multiple components." >> defense_systems_statistics.txt
    """
}