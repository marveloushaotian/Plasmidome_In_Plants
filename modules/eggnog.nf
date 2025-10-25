process EGGNOG_DIAMOND {
    tag "$sample_id"
    publishDir "${params.outdir}/eggnog_diamond", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(proteins_faa)
    val eggnog_db_path
    
    output:
    tuple val(sample_id), path("${sample_id}_diamond.emapper.seed_orthologs"), emit: diamond_results
    path "${sample_id}_diamond.log", emit: diamond_log
    
    script:
    """
    set -euo pipefail
    
    # Check if input file exists and has content
    if [ ! -s "${proteins_faa}" ]; then
        echo "[WARN] Input file ${proteins_faa} is empty or doesn't exist" >&2
        touch "${sample_id}_diamond.emapper.seed_orthologs"
        echo "No input file" > "${sample_id}_diamond.log"
        exit 0
    fi
    
    # Run eggNOG-mapper diamond search
    emapper.py -m diamond --no_annot --no_file_comments --cpu ${task.cpus} \\
        -i "${proteins_faa}" -o "${sample_id}_diamond" \\
        --data_dir "${eggnog_db_path}" -d none --itype proteins --override \\
        1>"${sample_id}_diamond.log" 2>&1 || {
        echo "[ERROR] eggNOG diamond search failed for ${sample_id}" >&2
        echo "Diamond search failed" > "${sample_id}_diamond.log"
        touch "${sample_id}_diamond.emapper.seed_orthologs"
        exit 0
    }
    """
}

process EGGNOG_ANNOTATION {
    tag "$sample_id"
    publishDir "${params.outdir}/eggnog_annotation", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(diamond_results)
    val eggnog_db_path
    
    output:
    tuple val(sample_id), path("${sample_id}_eggnog.emapper.annotations"), emit: annotations
    path "${sample_id}_annotation.log", emit: annotation_log
    
    script:
    """
    set -euo pipefail
    
    # Check if input file exists and has content
    if [ ! -s "${diamond_results}" ]; then
        echo "[WARN] Input file ${diamond_results} is empty or doesn't exist" >&2
        touch "${sample_id}_eggnog.emapper.annotations"
        echo "No input file" > "${sample_id}_annotation.log"
        exit 0
    fi
    
    # Run eggNOG-mapper annotation
    emapper.py --annotate_hits_table "${diamond_results}" \\
               --no_file_comments \\
               -o "${sample_id}_eggnog" \\
               --cpu ${task.cpus} \\
               --data_dir "${eggnog_db_path}" \\
               1>"${sample_id}_annotation.log" 2>&1 || {
        echo "[ERROR] eggNOG annotation failed for ${sample_id}" >&2
        echo "Annotation failed" > "${sample_id}_annotation.log"
        touch "${sample_id}_eggnog.emapper.annotations"
        exit 0
    }
    """
}

process EGGNOG_SUMMARY {
    tag "eggnog-summary"
    publishDir "${params.outdir}/eggnog_summary", mode: 'copy', overwrite: true
    
    input:
    path annotation_files
    
    output:
    path "eggnog_summary.tsv", emit: summary
    path "eggnog_statistics.txt", emit: statistics
    
    script:
    """
    set -euo pipefail
    
    echo "Sample_ID\\tTotal_Genes\\tAnnotated_Genes\\tAnnotation_Rate\\tCOG_Categories\\tKEGG_Pathways" > eggnog_summary.tsv
    echo "EggNOG Annotation Summary" > eggnog_statistics.txt
    echo "=========================" >> eggnog_statistics.txt
    echo "" >> eggnog_statistics.txt
    
    total_samples=0
    total_genes=0
    total_annotated=0
    
    for annotation_file in ${annotation_files}; do
        if [ -f "\$annotation_file" ] && [ -s "\$annotation_file" ]; then
            sample_id=\$(basename "\$annotation_file" _eggnog.emapper.annotations)
            
            # Count total genes (lines with data, excluding header)
            total_genes_count=\$(tail -n +2 "\$annotation_file" | wc -l)
            
            # Count annotated genes (non-empty COG_category or KEGG_ko)
            annotated_count=\$(tail -n +2 "\$annotation_file" | awk -F'\\t' 'NF > 5 && (\$6 != "" || \$7 != "")' | wc -l)
            
            # Calculate annotation rate
            if [ \$total_genes_count -gt 0 ]; then
                annotation_rate=\$(echo "scale=2; \$annotated_count * 100 / \$total_genes_count" | bc -l)
            else
                annotation_rate="0.00"
            fi
            
            # Count unique COG categories
            cog_categories=\$(tail -n +2 "\$annotation_file" | awk -F'\\t' '\$6 != "" {print \$6}' | sort | uniq | wc -l)
            
            # Count unique KEGG pathways
            kegg_pathways=\$(tail -n +2 "\$annotation_file" | awk -F'\\t' '\$7 != "" {print \$7}' | sort | uniq | wc -l)
            
            # Add to summary
            echo "\$sample_id\\t\$total_genes_count\\t\$annotated_count\\t\$annotation_rate\\t\$cog_categories\\t\$kegg_pathways" >> eggnog_summary.tsv
            
            # Add to statistics
            echo "Sample: \$sample_id" >> eggnog_statistics.txt
            echo "  Total genes: \$total_genes_count" >> eggnog_statistics.txt
            echo "  Annotated genes: \$annotated_count" >> eggnog_statistics.txt
            echo "  Annotation rate: \$annotation_rate%" >> eggnog_statistics.txt
            echo "  COG categories: \$cog_categories" >> eggnog_statistics.txt
            echo "  KEGG pathways: \$kegg_pathways" >> eggnog_statistics.txt
            echo "" >> eggnog_statistics.txt
            
            total_samples=\$((total_samples + 1))
            total_genes=\$((total_genes + total_genes_count))
            total_annotated=\$((total_annotated + annotated_count))
        fi
    done
    
    # Overall statistics
    if [ \$total_genes -gt 0 ]; then
        overall_rate=\$(echo "scale=2; \$total_annotated * 100 / \$total_genes" | bc -l)
    else
        overall_rate="0.00"
    fi
    
    echo "Overall Statistics:" >> eggnog_statistics.txt
    echo "  Total samples: \$total_samples" >> eggnog_statistics.txt
    echo "  Total genes: \$total_genes" >> eggnog_statistics.txt
    echo "  Total annotated: \$total_annotated" >> eggnog_statistics.txt
    echo "  Overall annotation rate: \$overall_rate%" >> eggnog_statistics.txt
    """
}