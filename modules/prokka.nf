process PROKKA_SANITIZE {
    tag "$sample_id"
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_sanitized.fasta"), emit: sanitized_fasta
    path "${sample_id}_contig_map.tsv", emit: contig_map
    
    script:
    """
    set -euo pipefail
    
    if [ ! -s "${fasta}" ]; then
        echo "[WARN] Input file ${fasta} is empty or doesn't exist" >&2
        touch "${sample_id}_sanitized.fasta"
        echo -e "original_id\\tnew_id" > "${sample_id}_contig_map.tsv"
        exit 0
    fi
    
    awk -v ID="${sample_id}" -v MAP="${sample_id}_contig_map.tsv" '
    BEGIN {
        FS="\\t"; OFS="\\t"; c=0;
        print "original_id","new_id" > MAP
    }
    /^>/{
        hdr = substr(\$0,2);               # remove ">"
        gsub(/[| ]+/, "_", hdr);          # replace spaces & pipes
        
        if (length(hdr) > 2000) { hdr = substr(hdr,1,2000) }
        
        c++
        new = sprintf("%s_c%05d", ID, c)
        print hdr, new >> MAP
        print ">" new
        next
    }
    { print }
    ' "${fasta}" > "${sample_id}_sanitized.fasta"
    
    if [ ! -s "${sample_id}_sanitized.fasta" ]; then
        echo "[ERROR] Sanitization failed for ${sample_id}" >&2
        exit 1
    fi
    """
}

process PROKKA_ANNOTATE {
    tag "$sample_id"
    publishDir "${params.outdir}/prokka/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'retry'
    maxRetries 2
    
    input:
    tuple val(sample_id), path(sanitized_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}.gff"), emit: gff optional true
    tuple val(sample_id), path("${sample_id}.faa"), emit: proteins optional true
    tuple val(sample_id), path("${sample_id}.ffn"), emit: genes optional true
    tuple val(sample_id), path("${sample_id}.gbk"), emit: genbank optional true
    path "${sample_id}.txt", emit: stats optional true
    path "*.log", emit: logs optional true
    
    script:
    """
    set -euo pipefail
    
    if [ ! -s "${sanitized_fasta}" ]; then
        echo "[WARN] Input file ${sanitized_fasta} is empty" >&2
        
        touch ${sample_id}.gff ${sample_id}.faa ${sample_id}.ffn ${sample_id}.gbk
        echo "Empty input for ${sample_id}" > ${sample_id}.txt
        echo "Empty input" > ${sample_id}.prokka.log
        exit 0
    fi
    
    seq_count=\$(grep -c "^>" "${sanitized_fasta}" || echo 0)
    if [ "\$seq_count" -eq 0 ]; then
        echo "[WARN] No sequences found in ${sanitized_fasta}" >&2
        touch ${sample_id}.gff ${sample_id}.faa ${sample_id}.ffn ${sample_id}.gbk
        echo "No sequences for ${sample_id}" > ${sample_id}.txt
        echo "No sequences" > ${sample_id}.prokka.log
        exit 0
    fi
    
    prokka \\
        --outdir . \\
        --prefix ${sample_id} \\
        --locustag ${sample_id} \\
        --kingdom ${params.prokka_kingdom} \\
        --cpus ${task.cpus} \\
        --force \\
        --quiet \\
        ${sanitized_fasta} \\
        2>&1 | tee ${sample_id}.prokka.log || {
        
        exitcode=\$?
        echo "[ERROR] Prokka failed for ${sample_id}, exit code: \$exitcode" >&2
        
        touch ${sample_id}.gff ${sample_id}.faa ${sample_id}.ffn ${sample_id}.gbk
        echo "Prokka failed for ${sample_id}" > ${sample_id}.txt
        exit 0
    }
    
    if [ ! -s "${sample_id}.gff" ] || [ ! -s "${sample_id}.faa" ]; then
        echo "[WARN] Prokka produced empty outputs for ${sample_id}" >&2
        
        touch ${sample_id}.gff ${sample_id}.faa ${sample_id}.ffn ${sample_id}.gbk
        echo "Incomplete annotation for ${sample_id}" > ${sample_id}.txt
    fi
    """
}

process PROKKA_SUMMARY {
    tag "prokka-summary"
    publishDir "${params.outdir}/prokka", mode: 'copy', overwrite: true
    
    input:
    path stat_files
    path protein_files
    
    output:
    path "annotation_summary.tsv", emit: summary
    path "all_proteins.faa", emit: combined_proteins
    path "proteins_for_eggnog.txt", emit: protein_list
    
    script:
    """
    set -euo pipefail
    
    echo -e "sample_id\\tcontigs\\tCDS\\trRNA\\ttRNA\\tmisc_RNA\\tstatus" > annotation_summary.tsv
    
    for stat in ${stat_files}; do
        if [ -f "\$stat" ]; then
            sample=\$(basename "\$stat" .txt)
            
            if grep -q "failed\\|Empty\\|No sequences" "\$stat" 2>/dev/null; then
                echo -e "\${sample}\\t0\\t0\\t0\\t0\\t0\\tfailed" >> annotation_summary.tsv
                continue
            fi
            
            contigs=\$(grep "contigs:" "\$stat" 2>/dev/null | awk '{print \$2}' || echo 0)
            cds=\$(grep "CDS:" "\$stat" 2>/dev/null | awk '{print \$2}' || echo 0)
            rrna=\$(grep "rRNA:" "\$stat" 2>/dev/null | awk '{print \$2}' || echo 0)
            trna=\$(grep "tRNA:" "\$stat" 2>/dev/null | awk '{print \$2}' || echo 0)
            misc=\$(grep "misc_RNA:" "\$stat" 2>/dev/null | awk '{print \$2}' || echo 0)
            
            echo -e "\${sample}\\t\${contigs:-0}\\t\${cds:-0}\\t\${rrna:-0}\\t\${trna:-0}\\t\${misc:-0}\\tsuccess" >> annotation_summary.tsv
        fi
    done
    
    > all_proteins.faa
    for faa in ${protein_files}; do
        if [ -s "\$faa" ]; then
            sample=\$(basename "\$faa" .faa)
            awk -v s="\$sample" '/^>/{print \$1"|"s; next} {print}' "\$faa" >> all_proteins.faa
        fi
    done
    
    ls *.faa 2>/dev/null | while read f; do
        [ -s "\$f" ] && echo "\$f"
    done | sort > proteins_for_eggnog.txt || touch proteins_for_eggnog.txt
    
    echo ""
    echo "Annotation Summary:"
    echo "=================="
    success=\$(grep -c "success\$" annotation_summary.tsv || echo 0)
    failed=\$(grep -c "failed\$" annotation_summary.tsv || echo 0)
    echo "Successfully annotated: \$success"
    echo "Failed: \$failed"
    """
}