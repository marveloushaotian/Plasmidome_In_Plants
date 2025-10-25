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
            grep "CDS" "\$gff" 2>/dev/null | awk -F'\\t' '{
                # Extract protein_id from attributes (column 9)
                match(\$9, /ID=([^;]+)/, id_arr)
                prot_id = id_arr[1]
                if (prot_id == "") {
                    match(\$9, /protein_id=([^;]+)/, pid_arr)
                    prot_id = pid_arr[1]
                }
                if (prot_id != "") {
                    print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$5"\\t"\$6"\\t"\$7"\\t"\$8"\\t"prot_id
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
    """
    #!/usr/bin/env python3
    ${file("${projectDir}/bin/merge_defense_systems.py").text}
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
        tail -n +2 ${defense_results} | awk -F'\\t' '\$3 != "" {print \$3}' | sort | uniq -c | sort -rn | head -20

        echo ""
        echo "PDC Systems:"
        echo "------------"
        pdc_count=\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$4 != "" {print \$4}' | wc -l)
        echo "Total PDC systems: \$pdc_count"

        echo ""
        echo "Antidefense Systems:"
        echo "--------------------"
        anti_count=\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$6 != "" {print \$6}' | wc -l)
        echo "Total antidefense systems: \$anti_count"

    } > defense_summary_report.txt

    # Create statistics table
    {
        echo -e "Category\\tCount"
        echo -e "Total_Systems\\t\$(tail -n +2 ${defense_results} | wc -l)"
        echo -e "Defense_Systems\\t\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$2 != "" {print \$2}' | wc -l)"
        echo -e "PDC_Systems\\t\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$4 != "" {print \$4}' | wc -l)"
        echo -e "Antidefense_Systems\\t\$(tail -n +2 ${defense_results} | awk -F'\\t' '\$6 != "" {print \$6}' | wc -l)"
    } > defense_statistics.tsv

    cat defense_summary_report.txt
    """
}
