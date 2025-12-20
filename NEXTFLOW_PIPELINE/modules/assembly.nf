process ASSEMBLY {
    tag "$sample_id"
    publishDir "${params.outdir}/assembly", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(reads)
    path adapters
    
    output:
    tuple val(sample_id), path("${sample_id}/contigs.fasta"), emit: contigs
    path "${sample_id}/scaffolds.fasta", emit: scaffolds
    path "${sample_id}/assembly_stats.txt", emit: stats
    
    script:
    """
    set -euo pipefail

    echo "[INFO] ASSEMBLY: sample=${sample_id}; cpus=${task.cpus}" >&2

    # BBDuk trimming
    bbduk.sh \\
        in1=${reads[0]} in2=${reads[1]} \\
        out1=${sample_id}_1.trimmed.fastq out2=${sample_id}_2.trimmed.fastq \\
        ref=${adapters} \\
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \\
        qtrim=r trimq=10 minlength=50 \\
        threads=${task.cpus}
    
    # Stage 1: Plasmid-focused assembly
    spades.py \\
        -1 ${sample_id}_1.trimmed.fastq -2 ${sample_id}_2.trimmed.fastq \\
        -o plasmid_assembly \\
        --only-assembler --plasmid \\
        -t ${task.cpus} || true
    
    # Stage 2: Final assembly
    if [ -f plasmid_assembly/scaffolds.fasta ]; then
        echo "[INFO] ${sample_id}: using plasmid_assembly/scaffolds.fasta as trusted-contigs" >&2
        spades.py \\
            -1 ${sample_id}_1.trimmed.fastq -2 ${sample_id}_2.trimmed.fastq \\
            -o ${sample_id} \\
            --only-assembler --isolate \\
            --trusted-contigs plasmid_assembly/scaffolds.fasta \\
            -t ${task.cpus}
    else
        echo "[WARN] ${sample_id}: no trusted-contigs found; running final assembly without them" >&2
        spades.py \\
            -1 ${sample_id}_1.trimmed.fastq -2 ${sample_id}_2.trimmed.fastq \\
            -o ${sample_id} \\
            --only-assembler --isolate \\
            -t ${task.cpus}
    fi
    
    # Generate statistics
    {
        echo "Sample: ${sample_id}"
        echo "Contigs: \$(grep -c '^>' ${sample_id}/contigs.fasta || echo 0)"
        echo "Total length: \$(awk '/^>/{next}{L+=length(\$0)}END{print (L?L:0)}' ${sample_id}/contigs.fasta)"
    } > ${sample_id}/assembly_stats.txt
    """
}
