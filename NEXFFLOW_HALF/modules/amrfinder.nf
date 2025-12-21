process AMRFINDER_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/amrfinder/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(protein_faa)

    output:
    tuple val(sample_id), path("${sample_id}_amrfinder_results.tsv"), emit: results
    tuple val(sample_id), path("${sample_id}_*"), emit: all_files, optional: true
    path "${sample_id}_amrfinder.log", emit: log

    script:
    """
    set -euo pipefail

    # Initialize log with basic info
    echo "[INFO] AMRFinder analysis for ${sample_id}" > "${sample_id}_amrfinder.log"
    echo "[INFO] Timestamp: \$(date)" >> "${sample_id}_amrfinder.log"
    echo "[INFO] Input file: ${protein_faa}" >> "${sample_id}_amrfinder.log"
    
    # Quick input validation
    if [ ! -f "${protein_faa}" ]; then
        echo "[ERROR] Input file does not exist: ${protein_faa}" >> "${sample_id}_amrfinder.log"
        echo "# AMRFinder failed - input file missing" > "${sample_id}_amrfinder_results.tsv"
        echo "Protein identifier\\tContig id\\tStart\\tStop\\tStrand\\tGene symbol\\tSequence name\\tScope\\tElement type\\tElement subtype\\tClass\\tSubclass\\tMethod\\tTarget length\\tReference sequence length\\t% Coverage of reference sequence\\t% Identity to reference sequence\\tAlignment length\\tAccession of closest sequence\\tName of closest sequence\\tHMM id\\tHMM description" >> "${sample_id}_amrfinder_results.tsv"
        exit 0
    fi
    
    file_size=\$(stat -c%s "${protein_faa}" 2>/dev/null || echo 0)
    seq_count=\$(grep -c "^>" "${protein_faa}" 2>/dev/null || echo 0)
    
    echo "[INFO] File size: \$file_size bytes" >> "${sample_id}_amrfinder.log"
    echo "[INFO] Protein sequences: \$seq_count" >> "${sample_id}_amrfinder.log"
    
    # Handle empty or invalid input
    if [ \$file_size -eq 0 ] || [ \$seq_count -eq 0 ]; then
        echo "[WARN] Empty or invalid input file" >> "${sample_id}_amrfinder.log"
        echo "# No proteins available for ${sample_id}" > "${sample_id}_amrfinder_results.tsv"
        echo "Protein identifier\\tContig id\\tStart\\tStop\\tStrand\\tGene symbol\\tSequence name\\tScope\\tElement type\\tElement subtype\\tClass\\tSubclass\\tMethod\\tTarget length\\tReference sequence length\\t% Coverage of reference sequence\\t% Identity to reference sequence\\tAlignment length\\tAccession of closest sequence\\tName of closest sequence\\tHMM id\\tHMM description" >> "${sample_id}_amrfinder_results.tsv"
        exit 0
    fi
    
    # Check AMRFinder availability
    if ! command -v amrfinder &> /dev/null; then
        echo "[ERROR] AMRFinder not found" >> "${sample_id}_amrfinder.log"
        echo "# AMRFinder not available" > "${sample_id}_amrfinder_results.tsv"
        echo "Protein identifier\\tContig id\\tStart\\tStop\\tStrand\\tGene symbol\\tSequence name\\tScope\\tElement type\\tElement subtype\\tClass\\tSubclass\\tMethod\\tTarget length\\tReference sequence length\\t% Coverage of reference sequence\\t% Identity to reference sequence\\tAlignment length\\tAccession of closest sequence\\tName of closest sequence\\tHMM id\\tHMM description" >> "${sample_id}_amrfinder_results.tsv"
        exit 0
    fi
    
    # Create output directory
    output_dir="amrfinder_\$\$"
    mkdir -p "\$output_dir"
    
    # Run AMRFinder with timeout
    echo "[INFO] Running AMRFinder..." >> "${sample_id}_amrfinder.log"
    timeout 1800 amrfinder -p "${protein_faa}" -o "\$output_dir/amrfinder_prot.out" > amrfinder.out 2> amrfinder.err
    exit_code=\$?
    
    echo "[INFO] AMRFinder exit code: \$exit_code" >> "${sample_id}_amrfinder.log"
    
    # Log outputs
    if [ -f amrfinder.out ]; then
        echo "[DEBUG] STDOUT:" >> "${sample_id}_amrfinder.log"
        cat amrfinder.out >> "${sample_id}_amrfinder.log"
    fi
    if [ -f amrfinder.err ]; then
        echo "[DEBUG] STDERR:" >> "${sample_id}_amrfinder.log"
        cat amrfinder.err >> "${sample_id}_amrfinder.log"
    fi
    
    # Process results
    if [ \$exit_code -eq 0 ] && [ -d "\$output_dir" ]; then
        echo "[INFO] AMRFinder completed successfully" >> "${sample_id}_amrfinder.log"
        
        # Copy output files
        for file in "\$output_dir"/*; do
            if [ -f "\$file" ]; then
                filename=\$(basename "\$file")
                cp "\$file" "${sample_id}_\$filename"
                echo "[INFO] Output: \$filename" >> "${sample_id}_amrfinder.log"
            fi
        done
        
        # Create main results file
        if [ -f "\$output_dir/amrfinder_prot.out" ]; then
            cp "\$output_dir/amrfinder_prot.out" "${sample_id}_amrfinder_results.tsv"
        else
            # AMRFinder succeeded but found no AMR genes
            echo "# No AMR genes found by AMRFinder" > "${sample_id}_amrfinder_results.tsv"
            echo "Protein identifier\\tContig id\\tStart\\tStop\\tStrand\\tGene symbol\\tSequence name\\tScope\\tElement type\\tElement subtype\\tClass\\tSubclass\\tMethod\\tTarget length\\tReference sequence length\\t% Coverage of reference sequence\\t% Identity to reference sequence\\tAlignment length\\tAccession of closest sequence\\tName of closest sequence\\tHMM id\\tHMM description" >> "${sample_id}_amrfinder_results.tsv"
        fi
    else
        echo "[WARN] AMRFinder failed or timed out" >> "${sample_id}_amrfinder.log"
        echo "# AMRFinder analysis failed" > "${sample_id}_amrfinder_results.tsv"
        echo "Protein identifier\\tContig id\\tStart\\tStop\\tStrand\\tGene symbol\\tSequence name\\tScope\\tElement type\\tElement subtype\\tClass\\tSubclass\\tMethod\\tTarget length\\tReference sequence length\\t% Coverage of reference sequence\\t% Identity to reference sequence\\tAlignment length\\tAccession of closest sequence\\tName of closest sequence\\tHMM id\\tHMM description" >> "${sample_id}_amrfinder_results.tsv"
    fi
    
    # Report final stats
    result_lines=\$(grep -v "^#" "${sample_id}_amrfinder_results.tsv" | tail -n +2 | wc -l 2>/dev/null || echo 0)
    echo "[INFO] Final result: \$result_lines AMR genes found" >> "${sample_id}_amrfinder.log"
    
    # Cleanup
    rm -rf "\$output_dir" amrfinder.out amrfinder.err
    """
}

process AMRFINDER_SUMMARY {
    tag "amrfinder-summary"
    publishDir "${params.outdir}/amrfinder_summary", mode: 'copy', overwrite: true

    input:
    path result_files

    output:
    path "amrfinder_summary.tsv", emit: summary
    path "amr_resistance_statistics.txt", emit: statistics
    
    when:
    result_files.size() > 0
    
    script:
    """
    set -euo pipefail
    
    # Create summary
    echo "Sample_ID\\tTotal_AMR\\tBeta-lactam\\tAminoglycoside\\tTetracycline\\tChloramphenicol\\tSulfonamide\\tQuinolone\\tMacrolide\\tOther\\tStatus" > amrfinder_summary.tsv
    
    # Statistics file
    echo "AMRFinder Resistance Gene Analysis Summary" > amr_resistance_statistics.txt
    echo "===========================================" >> amr_resistance_statistics.txt
    echo "Generated: \$(date)" >> amr_resistance_statistics.txt
    echo "" >> amr_resistance_statistics.txt
    
    total_samples=0
    successful_samples=0
    failed_samples=0
    samples_with_amr=0
    total_resistance_genes=0
    
    # AMR class counters
    beta_lactam_total=0
    aminoglycoside_total=0
    tetracycline_total=0
    chloramphenicol_total=0
    sulfonamide_total=0
    quinolone_total=0
    macrolide_total=0
    other_total=0
    
    for result_file in ${result_files}; do
        if [ -f "\$result_file" ]; then
            sample_id=\$(basename "\$result_file" _amrfinder_results.tsv)
            total_samples=\$((total_samples + 1))
            
            # Determine status
            if grep -q "^# AMRFinder failed\\|^# AMRFinder not available\\|^# No proteins" "\$result_file"; then
                status="Failed"
                failed_samples=\$((failed_samples + 1))
                amr_count=0
                beta_lactam=0
                aminoglycoside=0
                tetracycline=0
                chloramphenicol=0
                sulfonamide=0
                quinolone=0
                macrolide=0
                other=0
            else
                status="Success"
                successful_samples=\$((successful_samples + 1))
                amr_count=\$(grep -v "^#" "\$result_file" | tail -n +2 | wc -l 2>/dev/null || echo 0)
                
                if [ \$amr_count -gt 0 ]; then
                    samples_with_amr=\$((samples_with_amr + 1))
                    total_resistance_genes=\$((total_resistance_genes + amr_count))
                    
                    # Count resistance classes (based on common AMR gene patterns)
                    beta_lactam=\$(grep -v "^#" "\$result_file" | tail -n +2 | grep -i "beta-lactam\\|penicillin\\|cephalosporin\\|carbapenem\\|ampc\\|esbl" | wc -l 2>/dev/null || echo 0)
                    aminoglycoside=\$(grep -v "^#" "\$result_file" | tail -n +2 | grep -i "aminoglycoside\\|gentamicin\\|kanamycin\\|streptomycin\\|tobramycin" | wc -l 2>/dev/null || echo 0)
                    tetracycline=\$(grep -v "^#" "\$result_file" | tail -n +2 | grep -i "tetracycline\\|tet\\([A-Z]\\)" | wc -l 2>/dev/null || echo 0)
                    chloramphenicol=\$(grep -v "^#" "\$result_file" | tail -n +2 | grep -i "chloramphenicol\\|cat\\|cml\\|flo" | wc -l 2>/dev/null || echo 0)
                    sulfonamide=\$(grep -v "^#" "\$result_file" | tail -n +2 | grep -i "sulfonamide\\|sul\\|dhps\\|dhfr" | wc -l 2>/dev/null || echo 0)
                    quinolone=\$(grep -v "^#" "\$result_file" | tail -n +2 | grep -i "quinolone\\|fluoroquinolone\\|gyr\\|par\\|qnr" | wc -l 2>/dev/null || echo 0)
                    macrolide=\$(grep -v "^#" "\$result_file" | tail -n +2 | grep -i "macrolide\\|erythromycin\\|erm\\|mph\\|msr" | wc -l 2>/dev/null || echo 0)
                    
                    # Calculate other resistance genes
                    classified=\$((beta_lactam + aminoglycoside + tetracycline + chloramphenicol + sulfonamide + quinolone + macrolide))
                    if [ \$classified -le \$amr_count ]; then
                        other=\$((amr_count - classified))
                    else
                        other=0
                    fi
                    
                    # Add to totals
                    beta_lactam_total=\$((beta_lactam_total + beta_lactam))
                    aminoglycoside_total=\$((aminoglycoside_total + aminoglycoside))
                    tetracycline_total=\$((tetracycline_total + tetracycline))
                    chloramphenicol_total=\$((chloramphenicol_total + chloramphenicol))
                    sulfonamide_total=\$((sulfonamide_total + sulfonamide))
                    quinolone_total=\$((quinolone_total + quinolone))
                    macrolide_total=\$((macrolide_total + macrolide))
                    other_total=\$((other_total + other))
                else
                    beta_lactam=0
                    aminoglycoside=0
                    tetracycline=0
                    chloramphenicol=0
                    sulfonamide=0
                    quinolone=0
                    macrolide=0
                    other=0
                fi
            fi
            
            echo "\$sample_id\\t\$amr_count\\t\$beta_lactam\\t\$aminoglycoside\\t\$tetracycline\\t\$chloramphenicol\\t\$sulfonamide\\t\$quinolone\\t\$macrolide\\t\$other\\t\$status" >> amrfinder_summary.tsv
            
            echo "Sample \$sample_id: \$status (\$amr_count resistance genes)" >> amr_resistance_statistics.txt
        fi
    done
    
    echo "" >> amr_resistance_statistics.txt
    echo "Overall Statistics:" >> amr_resistance_statistics.txt
    echo "- Total samples: \$total_samples" >> amr_resistance_statistics.txt
    echo "- Successful: \$successful_samples" >> amr_resistance_statistics.txt
    echo "- Failed: \$failed_samples" >> amr_resistance_statistics.txt
    echo "- With AMR genes: \$samples_with_amr" >> amr_resistance_statistics.txt
    echo "- Total resistance genes: \$total_resistance_genes" >> amr_resistance_statistics.txt
    echo "" >> amr_resistance_statistics.txt
    echo "Resistance Class Breakdown:" >> amr_resistance_statistics.txt
    echo "- Beta-lactam: \$beta_lactam_total" >> amr_resistance_statistics.txt
    echo "- Aminoglycoside: \$aminoglycoside_total" >> amr_resistance_statistics.txt
    echo "- Tetracycline: \$tetracycline_total" >> amr_resistance_statistics.txt
    echo "- Chloramphenicol: \$chloramphenicol_total" >> amr_resistance_statistics.txt
    echo "- Sulfonamide: \$sulfonamide_total" >> amr_resistance_statistics.txt
    echo "- Quinolone: \$quinolone_total" >> amr_resistance_statistics.txt
    echo "- Macrolide: \$macrolide_total" >> amr_resistance_statistics.txt
    echo "- Other: \$other_total" >> amr_resistance_statistics.txt
    
    if [ \$successful_samples -gt 0 ]; then
        avg=\$(echo "\$total_resistance_genes \$successful_samples" | awk '{printf "%.2f", \$1 / \$2}')
        echo "- Average per successful sample: \$avg" >> amr_resistance_statistics.txt
    fi
    """
}
