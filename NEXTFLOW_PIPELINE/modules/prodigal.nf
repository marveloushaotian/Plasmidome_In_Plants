process PRODIGAL_PREDICT {
    tag "$sample_id"
    publishDir "${params.outdir}/prodigal", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}.gff3"), emit: gff
    tuple val(sample_id), path("${sample_id}.fna"), emit: genes
    tuple val(sample_id), path("${sample_id}.faa"), emit: proteins
    path "${sample_id}.log", emit: log
    
    script:
    """
    set -euo pipefail
    
    # Check if input file exists and has content
    if [ ! -s "${fasta}" ]; then
        echo "[WARN] Input file ${fasta} is empty or doesn't exist" >&2
        touch "${sample_id}.gff3"
        touch "${sample_id}.fna"
        touch "${sample_id}.faa"
        echo "No input file" > "${sample_id}.log"
        exit 0
    fi
    
    # Run Prodigal gene prediction
    prodigal -i "${fasta}" -f gff -o "${sample_id}.gff3" -d "${sample_id}.fna" -a "${sample_id}.faa" &> "${sample_id}.log" || {
        echo "[ERROR] Prodigal failed for ${sample_id}" >&2
        echo "Prodigal prediction failed" > "${sample_id}.log"
        touch "${sample_id}.gff3"
        touch "${sample_id}.fna"
        touch "${sample_id}.faa"
        exit 0
    }
    
    echo "[INFO] Prodigal completed for ${sample_id}" >&2
    """
}

process PRODIGAL_SUMMARY {
    tag "prodigal-summary"
    publishDir "${params.outdir}/prodigal_summary", mode: 'copy', overwrite: true
    
    input:
    path gff_files
    path gene_files
    path protein_files
    
    output:
    path "prodigal_summary.tsv", emit: summary
    path "prodigal_statistics.txt", emit: statistics
    
    script:
    """
    set -euo pipefail
    
    echo "Sample_ID\\tTotal_Genes\\tTotal_Proteins\\tAverage_Gene_Length\\tGC_Content" > prodigal_summary.tsv
    echo "Prodigal Gene Prediction Summary" > prodigal_statistics.txt
    echo "=================================" >> prodigal_statistics.txt
    echo "" >> prodigal_statistics.txt
    
    total_samples=0
    total_genes=0
    total_proteins=0
    
    for gff_file in ${gff_files}; do
        if [ -f "\$gff_file" ] && [ -s "\$gff_file" ]; then
            sample_id=\$(basename "\$gff_file" .gff3)
            
            # Count genes from GFF file
            gene_count=\$(grep -c "^[^#]" "\$gff_file" 2>/dev/null || echo 0)
            
            # Count proteins from corresponding protein file
            protein_file="\${gff_file%.gff3}.faa"
            protein_count=0
            if [ -f "\$protein_file" ] && [ -s "\$protein_file" ]; then
                protein_count=\$(grep -c "^>" "\$protein_file" 2>/dev/null || echo 0)
            fi
            
            # Calculate average gene length from GFF
            avg_length=0
            if [ \$gene_count -gt 0 ]; then
                avg_length=\$(awk -F'\\t' 'BEGIN{sum=0; count=0} /^[^#]/{sum+=\$5-\$4+1; count++} END{if(count>0) printf "%.1f", sum/count; else print "0"}' "\$gff_file" 2>/dev/null || echo "0")
            fi
            
            # Get GC content from log file
            log_file="\${gff_file%.gff3}.log"
            gc_content="0.0"
            if [ -f "\$log_file" ]; then
                gc_content=\$(grep "GC content:" "\$log_file" | awk '{print \$3}' | sed 's/%//' 2>/dev/null || echo "0.0")
            fi
            
            # Add to summary
            echo "\$sample_id\\t\$gene_count\\t\$protein_count\\t\$avg_length\\t\$gc_content" >> prodigal_summary.tsv
            
            # Add to statistics
            echo "Sample: \$sample_id" >> prodigal_statistics.txt
            echo "  Total genes: \$gene_count" >> prodigal_statistics.txt
            echo "  Total proteins: \$protein_count" >> prodigal_statistics.txt
            echo "  Average gene length: \$avg_length bp" >> prodigal_statistics.txt
            echo "  GC content: \$gc_content%" >> prodigal_statistics.txt
            echo "" >> prodigal_statistics.txt
            
            total_samples=\$((total_samples + 1))
            total_genes=\$((total_genes + gene_count))
            total_proteins=\$((total_proteins + protein_count))
        fi
    done
    
    # Overall statistics
    if [ \$total_genes -gt 0 ]; then
        avg_genes_per_sample=\$(echo "scale=2; \$total_genes / \$total_samples" | bc -l)
    else
        avg_genes_per_sample="0.00"
    fi
    
    echo "Overall Statistics:" >> prodigal_statistics.txt
    echo "  Total samples: \$total_samples" >> prodigal_statistics.txt
    echo "  Total genes: \$total_genes" >> prodigal_statistics.txt
    echo "  Total proteins: \$total_proteins" >> prodigal_statistics.txt
    echo "  Average genes per sample: \$avg_genes_per_sample" >> prodigal_statistics.txt
    """
}