process MOBSUITE_RECON {
    tag "$sample_id"
    publishDir "${params.outdir}/mobsuite_recon/${sample_id}", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(fasta)
    
    output:
    tuple val(sample_id), path("*.fasta"), emit: plasmid_fastas optional true
    tuple val(sample_id), path("${sample_id}_contig_report.txt"), emit: contig_report optional true
    path "*.txt", emit: all_reports optional true
    
    script:
    """
    set -euo pipefail
    
    if [ ! -s "${fasta}" ]; then
        echo "[WARN] Input file ${fasta} is empty or doesn't exist" >&2
        echo "sample_id\tstatus" > ${sample_id}_contig_report.txt
        echo "${sample_id}\tempty_input" >> ${sample_id}_contig_report.txt
        exit 0
    fi
    
    mob_recon \\
        --infile ${fasta} \\
        --outdir . \\
        --force \\
        --num_threads ${task.cpus} 2>&1 | tee mob_recon.log || {
        echo "[ERROR] mob_recon failed for ${sample_id}, exit code: \$?" >&2
        
        echo "sample_id\tstatus" > ${sample_id}_contig_report.txt
        echo "${sample_id}\tfailed" >> ${sample_id}_contig_report.txt
        exit 0
    }
    
    if [ -f "contig_report.txt" ]; then
        mv contig_report.txt ${sample_id}_contig_report.txt
    else
        echo "[WARN] No contig_report.txt generated for ${sample_id}" >&2
        echo "sample_id\tstatus" > ${sample_id}_contig_report.txt
        echo "${sample_id}\tno_output" >> ${sample_id}_contig_report.txt
    fi
    
    if [ -f "chromosome.fasta" ]; then
        mv chromosome.fasta ${sample_id}_chromosome.fasta 2>/dev/null || true
    fi
    
    for plasmid in plasmid_*.fasta; do
        if [ -f "\$plasmid" ]; then
            base=\$(basename "\$plasmid" .fasta)
            mv "\$plasmid" ${sample_id}_\${base}.fasta 2>/dev/null || true
        fi
    done
    """
}

process MOBSUITE_TYPER {
    tag "$sample_id"
    publishDir "${params.outdir}/mobsuite_typer", mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_mobtyper.tsv"), emit: typing_result optional true
    
    script:
    """
    set -euo pipefail
    
    if [ ! -s "${fasta}" ]; then
        echo "[WARN] Input file ${fasta} is empty or doesn't exist" >&2
        echo -e "file_id\tnum_contigs\ttotal_length\tmobility" > ${sample_id}_mobtyper.tsv
        echo -e "${sample_id}\t0\t0\tNA" >> ${sample_id}_mobtyper.tsv
        exit 0
    fi
    
    mob_typer \\
        --multi \\
        --infile ${fasta} \\
        --out_file ${sample_id}_mobtyper.tsv 2>&1 | tee mob_typer.log || {
        echo "[ERROR] mob_typer failed for ${sample_id}, exit code: \$?" >&2
        
        echo -e "file_id\tnum_contigs\ttotal_length\tmobility" > ${sample_id}_mobtyper.tsv
        echo -e "${sample_id}\t0\t0\tfailed" >> ${sample_id}_mobtyper.tsv
        exit 0
    }
    
    if [ ! -f "${sample_id}_mobtyper.tsv" ]; then
        echo "[WARN] No output generated for ${sample_id}" >&2
        echo -e "file_id\tnum_contigs\ttotal_length\tmobility" > ${sample_id}_mobtyper.tsv
        echo -e "${sample_id}\t0\t0\tno_output" >> ${sample_id}_mobtyper.tsv
    fi
    """
}

process MOBSUITE_SUMMARY {
    tag "mobsuite-summary"
    publishDir "${params.outdir}/mobsuite", mode: 'copy', overwrite: true
    
    input:
    path typing_files
    path contig_reports
    
    output:
    path "mobsuite_summary.tsv", emit: summary
    path "plasmid_statistics.txt", emit: stats
    
    script:
    """
    set -euo pipefail
    
    echo "Typing files received:" >&2
    ls -la *_mobtyper.tsv 2>/dev/null || echo "No typing files" >&2
    
    echo "Contig reports received:" >&2
    ls -la *_contig_report.txt 2>/dev/null || echo "No contig reports" >&2
    
    echo -e "sample_id\\tcontig_id\\tmobility\\tinc_groups\\tmash_hit" > mobsuite_summary.tsv
    
    for tsv in *_mobtyper.tsv; do
        if [ -s "\$tsv" ]; then
            sample=\$(basename "\$tsv" _mobtyper.tsv)
            tail -n +2 "\$tsv" 2>/dev/null | grep -v "^failed" | grep -v "^no_output" | awk -v s="\$sample" 'NF>0 {print s"\\t"\$0}' >> mobsuite_summary.tsv || true
        fi
    done
    
    {
        echo "MOB-suite Analysis Summary"
        echo "=========================="
        echo "Total samples processed: \$(ls *_mobtyper.tsv 2>/dev/null | wc -l || echo 0)"
        
        success_count=\$(tail -n +2 mobsuite_summary.tsv | wc -l || echo 0)
        echo "Total plasmids/contigs found: \$success_count"
        
        samples_with_plasmids=\$(tail -n +2 mobsuite_summary.tsv | cut -f1 | sort -u | wc -l || echo 0)
        echo "Samples with mobile elements: \$samples_with_plasmids"
        
        echo ""
        echo "Plasmid mobility distribution:"
        if [ -s mobsuite_summary.tsv ]; then
            tail -n +2 mobsuite_summary.tsv | awk '{mob[\$3]++} END {for(m in mob) print m": "mob[m]}' | sort || echo "No mobile elements detected"
        else
            echo "No mobile elements detected"
        fi
        
        echo ""
        echo "Samples without mobile elements or failed:"
        for report in *_contig_report.txt; do
            if [ -s "\$report" ]; then
                sample=\$(basename "\$report" _contig_report.txt)
                if ! grep -q "^\${sample}\\t" mobsuite_summary.tsv 2>/dev/null; then
                    echo "\$sample"
                fi
            fi
        done | sort -u | head -20
    } > plasmid_statistics.txt
    
    echo "[INFO] Summary complete. Check plasmid_statistics.txt for details." >&2
    """
}
