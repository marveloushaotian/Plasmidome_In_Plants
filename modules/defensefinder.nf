process DEFENSEFINDER_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/defensefinder/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(protein_faa)
    
    output:
    tuple val(sample_id), path("${sample_id}_defensefinder.csv"), emit: results
    tuple val(sample_id), path("${sample_id}_*"), emit: all_files
    path "${sample_id}_defensefinder.log", emit: log
    
    script:
    """
    set -euo pipefail
    
    echo "[INFO] Starting DefenseFinder annotation for ${sample_id}" > "${sample_id}_defensefinder.log"
    echo "[INFO] Timestamp: \$(date)" >> "${sample_id}_defensefinder.log"
    echo "[INFO] Working directory: \$(pwd)" >> "${sample_id}_defensefinder.log"
    echo "[INFO] Using protein sequence mode" >> "${sample_id}_defensefinder.log"
    
    if [ ! -s "${protein_faa}" ]; then
        echo "[WARN] Input protein file ${protein_faa} is empty or doesn't exist" >&2
        echo "[WARN] Input protein file ${protein_faa} is empty or doesn't exist" >> "${sample_id}_defensefinder.log"
        touch "${sample_id}_defensefinder.csv"
        exit 0
    fi
    
    if ! command -v defense-finder &> /dev/null; then
        echo "[ERROR] defense-finder not found in PATH" >&2
        echo "[ERROR] defense-finder not found in PATH" >> "${sample_id}_defensefinder.log"
        touch "${sample_id}_defensefinder.csv"
        exit 0
    fi
    
    echo "[INFO] Protein file: ${protein_faa}" >> "${sample_id}_defensefinder.log"
    echo "[INFO] Sample ID: ${sample_id}" >> "${sample_id}_defensefinder.log"
    
    echo "[INFO] DefenseFinder version:" >> "${sample_id}_defensefinder.log"
    defense-finder --version >> "${sample_id}_defensefinder.log" 2>&1 || echo "[WARN] Version check failed" >> "${sample_id}_defensefinder.log"
    
    echo "[INFO] Input file information:" >> "${sample_id}_defensefinder.log"
    echo "  File size: \$(stat -c%s '${protein_faa}' 2>/dev/null || echo 'unknown') bytes" >> "${sample_id}_defensefinder.log"
    echo "  Protein sequences count: \$(grep -c '^>' '${protein_faa}' 2>/dev/null || echo 'unknown')" >> "${sample_id}_defensefinder.log"
    
    mkdir -p defensefinder_output
    echo "[INFO] Created output directory: defensefinder_output" >> "${sample_id}_defensefinder.log"
    
    echo "[INFO] Running DefenseFinder command:" >> "${sample_id}_defensefinder.log"
    echo "defense-finder run -a --db-type gembase ${protein_faa} -o defensefinder_output" >> "${sample_id}_defensefinder.log"
    echo "[INFO] Starting DefenseFinder execution at: \$(date)" >> "${sample_id}_defensefinder.log"
    
    defense-finder run -a --db-type gembase "${protein_faa}" -o defensefinder_output > defensefinder_stdout.log 2> defensefinder_stderr.log
    exit_code=\$?
    
    echo "[INFO] DefenseFinder finished at: \$(date)" >> "${sample_id}_defensefinder.log"
    echo "[INFO] DefenseFinder exit code: \$exit_code" >> "${sample_id}_defensefinder.log"
    
    echo "[INFO] DefenseFinder stdout:" >> "${sample_id}_defensefinder.log"
    cat defensefinder_stdout.log >> "${sample_id}_defensefinder.log" 2>/dev/null || echo "No stdout captured" >> "${sample_id}_defensefinder.log"
    echo "[INFO] DefenseFinder stderr:" >> "${sample_id}_defensefinder.log"
    cat defensefinder_stderr.log >> "${sample_id}_defensefinder.log" 2>/dev/null || echo "No stderr captured" >> "${sample_id}_defensefinder.log"
    
    echo "[INFO] Checking for output files:" >> "${sample_id}_defensefinder.log"
    ls -la defensefinder_output/ >> "${sample_id}_defensefinder.log" 2>&1 || echo "No defensefinder_output directory found" >> "${sample_id}_defensefinder.log"
    
    tsv_found=false
    tsv_source=""
    
    if [ -f "defensefinder_output/${sample_id}_defense_finder_systems.tsv" ] && [ -s "defensefinder_output/${sample_id}_defense_finder_systems.tsv" ]; then
        tsv_source="defensefinder_output/${sample_id}_defense_finder_systems.tsv"
        tsv_found=true
    else
        tsv_file=\$(find defensefinder_output/ -name "*_defense_finder_systems.tsv" -type f -size +0c | head -1)
        if [ -n "\$tsv_file" ] && [ -f "\$tsv_file" ]; then
            tsv_source="\$tsv_file"
            tsv_found=true
        fi
    fi
    
    if [ "\$tsv_found" = true ] && [ \$exit_code -eq 0 ]; then
        cp "\$tsv_source" "${sample_id}_defensefinder.csv"
        
        output_size=\$(stat -c%s "${sample_id}_defensefinder.csv" 2>/dev/null || echo 0)
        systems_count=\$(tail -n +2 "${sample_id}_defensefinder.csv" | wc -l 2>/dev/null || echo 0)
        
        echo "[INFO] DefenseFinder completed successfully for ${sample_id}" >&2
        echo "[INFO] Output file size: \$output_size bytes" >> "${sample_id}_defensefinder.log"
        echo "[INFO] Defense systems found: \$systems_count" >> "${sample_id}_defensefinder.log"
        
        for file in defensefinder_output/*; do
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
        echo "[WARN] DefenseFinder failed or no output found for ${sample_id}" >&2
        echo "[WARN] DefenseFinder exit code: \$exit_code" >> "${sample_id}_defensefinder.log"
        touch "${sample_id}_defensefinder.csv"
    fi
    
    rm -f defensefinder_stdout.log defensefinder_stderr.log
    
    echo "[INFO] Process completed at: \$(date)" >> "${sample_id}_defensefinder.log"
    """
}

process DEFENSEFINDER_SUMMARY {
    tag "defensefinder-summary"
    publishDir "${params.outdir}/defensefinder_summary", mode: 'copy', overwrite: true
    
    input:
    path result_files
    
    output:
    path "defensefinder_summary.tsv", emit: summary
    path "defense_systems_statistics.txt", emit: statistics
    
    script:
    """
    set -euo pipefail
    
    echo "Sample_ID\\tTotal_Systems\\tCRISPR\\tRestriction_Modification\\tAbi\\tToxin_Antitoxin\\tOther\\tFile_Size_Bytes" > defensefinder_summary.tsv
    
    echo "DefenseFinder Defense Systems Summary" > defense_systems_statistics.txt
    echo "====================================" >> defense_systems_statistics.txt
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
            sample_id=\$(basename "\$result_file" _defensefinder.csv)
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
                
                crispr_count=\$(tail -n +2 "\$result_file" | cut -d'\t' -f3 | grep -i "crispr\\|cas\\|type.*i\\|type.*ii\\|type.*iii\\|type.*iv\\|type.*v\\|type.*vi" | wc -l 2>/dev/null || echo 0)
                
                rm_count=\$(tail -n +2 "\$result_file" | cut -d'\t' -f3 | grep -i "restriction\\|modification\\|rm\\|hsdm\\|hsdr\\|hsds" | wc -l 2>/dev/null || echo 0)
                
                abi_count=\$(tail -n +2 "\$result_file" | cut -d'\t' -f3 | grep -i "abortive\\|abi\\|abort" | wc -l 2>/dev/null || echo 0)
                
                ta_count=\$(tail -n +2 "\$result_file" | cut -d'\t' -f3 | grep -i "toxin\\|antitoxin\\|mazf\\|rele\\|hipa\\|vapb\\|vapc" | wc -l 2>/dev/null || echo 0)
                
                classified_count=\$((crispr_count + rm_count + abi_count + ta_count))
                if [ \$classified_count -le \$total_systems_count ]; then
                    other_count=\$((total_systems_count - classified_count))
                else
                    other_count=0
                fi
            fi
            
            echo "\$sample_id\\t\$total_systems_count\\t\$crispr_count\\t\$rm_count\\t\$abi_count\\t\$ta_count\\t\$other_count\\t\$file_size" >> defensefinder_summary.tsv
            
            echo "Sample: \$sample_id" >> defense_systems_statistics.txt
            echo "  Total defense systems: \$total_systems_count" >> defense_systems_statistics.txt
            if [ \$total_systems_count -gt 0 ]; then
                echo "  System breakdown:" >> defense_systems_statistics.txt
                tail -n +2 "\$result_file" | cut -d'\t' -f3 | sort | uniq -c | sed 's/^/    /' >> defense_systems_statistics.txt
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
    echo "Note: DefenseFinder analysis provides system-level identification." >> defense_systems_statistics.txt
    echo "Individual systems may contain multiple components." >> defense_systems_statistics.txt
    """
}
