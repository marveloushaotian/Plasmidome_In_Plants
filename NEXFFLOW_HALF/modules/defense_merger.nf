process COLLECT_GFF_FILES {
    tag "collect-gff"
    publishDir "${params.outdir}/defense_merger/intermediate", mode: 'copy', overwrite: true

    input:
    path gff_files

    output:
    path "all_transformed_gff.tsv", emit: merged_gff

    script:
    """
    set -euo pipefail

    # Create header
    echo -e "seqid\\tsource\\ttype\\tstart\\tend\\tscore\\tstrand\\tphase\\tprot_id" > all_transformed_gff.tsv

    # Combine all GFF files
    for gff in ${gff_files}; do
        if [ -s "\$gff" ]; then
            # Extract CDS features and convert to TSV format
            # Prodigal protein IDs: "1_31" means contig 1, gene 31
            # Replace contig number with full seqid: "1_31" -> "NODE_1_..._31"
            grep "CDS" "\$gff" 2>/dev/null | awk -F'\\t' '{
                # Extract protein_id from attributes (column 9)
                match(\$9, /ID=([^;]+)/, id_arr)
                prot_id = id_arr[1]
                if (prot_id == "") {
                    match(\$9, /protein_id=([^;]+)/, pid_arr)
                    prot_id = pid_arr[1]
                }
                if (prot_id != "") {
                    # Remove contig number from prot_id (everything before first underscore)
                    # e.g., "1_31" -> "31"
                    sub(/^[^_]+_/, "", prot_id)
                    # Create full protein ID: seqid + "_" + gene_number
                    full_prot_id = \$1 "_" prot_id
                    print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"\$8"\\t"full_prot_id
                }
            }' >> all_transformed_gff.tsv || true
        fi
    done

    # Check if we got any data
    line_count=\$(wc -l < all_transformed_gff.tsv)
    if [ "\$line_count" -le 1 ]; then
        echo "[WARN] No CDS features found in GFF files" >&2
    fi
    """
}

process COLLECT_PADLOC_RESULTS {
    tag "collect-padloc"
    publishDir "${params.outdir}/defense_merger/intermediate", mode: 'copy', overwrite: true

    input:
    path padloc_files

    output:
    path "final_padloc.tsv", emit: merged_padloc

    script:
    """
    set -euo pipefail

    # Combine all PADLOC results
    first_file=true
    > final_padloc.tsv

    for file in ${padloc_files}; do
        if [ -s "\$file" ]; then
            if [ "\$first_file" = true ]; then
                cat "\$file" >> final_padloc.tsv
                first_file=false
            else
                # Skip header for subsequent files
                tail -n +2 "\$file" >> final_padloc.tsv || true
            fi
        fi
    done

    # Create minimal file if empty
    if [ ! -s final_padloc.tsv ]; then
        echo -e "seqid\\tsystem\\tsystem.number\\ttarget.name" > final_padloc.tsv
    fi
    """
}

process COLLECT_DEFENSEFINDER_RESULTS {
    tag "collect-defensefinder"
    publishDir "${params.outdir}/defense_merger/intermediate", mode: 'copy', overwrite: true

    input:
    path defensefinder_files

    output:
    path "final_defensefinder.tsv", emit: merged_defensefinder

    script:
    """
    set -euo pipefail

    # Combine all DefenseFinder results
    first_file=true
    > final_defensefinder.tsv

    for file in ${defensefinder_files}; do
        if [ -s "\$file" ]; then
            if [ "\$first_file" = true ]; then
                cat "\$file" >> final_defensefinder.tsv
                first_file=false
            else
                tail -n +2 "\$file" >> final_defensefinder.tsv || true
            fi
        fi
    done

    # Create minimal file if empty
    if [ ! -s final_defensefinder.tsv ]; then
        echo -e "sys_id\\ttype\\tsubtype\\tprotein_in_syst\\tname_of_profiles_in_sys\\tactivity\\tsys_beg" > final_defensefinder.tsv
    fi
    """
}

process COLLECT_CCTYPER_RESULTS {
    tag "collect-cctyper"
    publishDir "${params.outdir}/defense_merger/intermediate", mode: 'copy', overwrite: true

    input:
    path cctyper_files

    output:
    path "final_cctyper.tsv", emit: merged_cctyper

    script:
    """
    set -euo pipefail

    # Combine all CCTyper results
    first_file=true
    > final_cctyper.tsv

    for file in ${cctyper_files}; do
        if [ -s "\$file" ]; then
            if [ "\$first_file" = true ]; then
                cat "\$file" >> final_cctyper.tsv
                first_file=false
            else
                tail -n +2 "\$file" >> final_cctyper.tsv || true
            fi
        fi
    done

    # Create minimal file if empty
    if [ ! -s final_cctyper.tsv ]; then
        echo -e "Prediction\\tGenes\\tProt_IDs\\tOperon" > final_cctyper.tsv
    fi
    """
}

process MERGE_DEFENSE_SYSTEMS {
    tag "merge-defense-systems"
    publishDir "${params.outdir}/defense_merger", mode: 'copy', overwrite: true

    input:
    path merged_gff
    path merged_padloc
    path merged_defensefinder
    path merged_cctyper
    path mapping_file

    output:
    path "all.merged", emit: all_merged
    path "all.merged.unified", emit: unified_merged
    path "defense_results.tsv", emit: defense_summary
    path "merge_stats.txt", emit: stats

    script:
    def merge_script = file("${projectDir}/bin/merge_defense_systems.py")
    def merge_signature = merge_script.text.hashCode()
    """
    set -euo pipefail

    echo "merge_defense_systems.py signature: ${merge_signature}"

    # Run merge_defense_systems.py with proper arguments
    python3 ${merge_script} \
        --mapping ${mapping_file} \
        --gff ${merged_gff} \
        --padloc ${merged_padloc} \
        --defensefinder ${merged_defensefinder} \
        --cctyper ${merged_cctyper}
    
    # Generate merge statistics
    echo "Merge statistics:" > merge_stats.txt
    echo "Input files:" >> merge_stats.txt
    echo "  GFF: ${merged_gff}" >> merge_stats.txt
    echo "  PADLOC: ${merged_padloc}" >> merge_stats.txt
    echo "  DefenseFinder: ${merged_defensefinder}" >> merge_stats.txt
    echo "  CCTyper: ${merged_cctyper}" >> merge_stats.txt
    echo "  Mapping: ${mapping_file}" >> merge_stats.txt
    echo "" >> merge_stats.txt
    echo "Output files:" >> merge_stats.txt
    if [ -f "all.merged" ]; then
        echo "  all.merged: \$(wc -l < all.merged) lines" >> merge_stats.txt
    fi
    if [ -f "all.merged.unified" ]; then
        echo "  all.merged.unified: \$(wc -l < all.merged.unified) lines" >> merge_stats.txt
    fi
    if [ -f "defense_results.tsv" ]; then
        echo "  defense_results.tsv: \$(wc -l < defense_results.tsv) lines" >> merge_stats.txt
    fi
    """
}

process DEFENSE_SUMMARY {
    tag "defense-summary"
    publishDir "${params.outdir}/defense_merger", mode: 'copy', overwrite: true

    input:
    path defense_results

    output:
    path "defense_summary_report.txt", emit: report
    path "defense_statistics.tsv", emit: statistics

    script:
    """
    set -euo pipefail

    {
        echo "Defense Systems Analysis Summary"
        echo "================================="
        echo ""
        echo "Total defense systems identified: \$(tail -n +2 ${defense_results} | wc -l)"
        echo ""

        echo "By Type:"
        echo "--------"
        tail -n +2 ${defense_results} | awk -F'\\t' '\$2 != "" {print \$2}' | sort | uniq -c | sort -rn

        echo ""
        echo "By Subtype (top 20):"
        echo "--------------------"
        # 临时禁用pipefail以避免head命令的SIGPIPE问题
        set +o pipefail
        tail -n +2 ${defense_results} | awk -F'\\t' '\$3 != "" {print \$3}' | sort | uniq -c | sort -rn | head -20
        set -o pipefail

        echo ""
        echo "PDC Systems:"
        echo "------------"
        pdc_count=\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$4 != "" {print \$4}' | wc -l)
        echo "Total PDC systems: \$pdc_count"

        echo ""
        echo "Antidefense Systems:"
        echo "--------------------"
        anti_count=\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$18 != "" {print \$18}' | wc -l)
        echo "Total antidefense systems: \$anti_count"

    } > defense_summary_report.txt

    # Create statistics table
    {
        echo -e "Category\\tCount"
        echo -e "Total_Systems\\t\$(tail -n +2 ${defense_results} | wc -l)"
        echo -e "Defense_Systems\\t\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$2 != "" {print \$2}' | wc -l)"
        echo -e "PDC_Systems\\t\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$4 != "" {print \$4}' | wc -l)"
        echo -e "Antidefense_Systems\\t\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$18 != "" {print \$18}' | wc -l)"
    } > defense_statistics.tsv

    cat defense_summary_report.txt
    """
}
