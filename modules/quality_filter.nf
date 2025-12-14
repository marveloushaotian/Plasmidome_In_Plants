process QUALITY_FILTER {
    tag "$sample_id"
    publishDir "${params.outdir}/filtered", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(fasta), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered.fasta"), emit: filtered_fasta
    path "${sample_id}_coverage_stats.txt", emit: coverage_stats
    
    script:
    """
    set -euo pipefail

    echo "[INFO] QUALITY_FILTER: sample=${sample_id}; cpus=${task.cpus}; len>=${params.length_threshold}; cov>=${params.coverage_threshold}" >&2

    reformat.sh \\
        in=${fasta} \\
        out=${sample_id}_length_filtered.fasta \\
        minlength=${params.length_threshold} \\
        ow=t
    
    bbmap.sh \\
        ref=${sample_id}_length_filtered.fasta \\
        in1=${reads[0]} in2=${reads[1]} \\
        out=${sample_id}_mapped.sam \\
        threads=${task.cpus} \\
        nodisk
    
    pileup.sh \\
        in=${sample_id}_mapped.sam \\
        covstats=${sample_id}_contig_covstats.txt
    
    awk -v threshold=${params.coverage_threshold} \\
        'NR>1 && \$2 >= threshold {print \$1}' ${sample_id}_contig_covstats.txt > keep_ids.txt
    
    filterbyname.sh \\
        in=${sample_id}_length_filtered.fasta \\
        out=${sample_id}_filtered.fasta \\
        include=t \\
        names=keep_ids.txt \\
        overwrite=t
    
    cp ${sample_id}_contig_covstats.txt ${sample_id}_coverage_stats.txt
    """
}
