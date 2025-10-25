process GTDBTK {
    tag "gtdbtk-batch"
    publishDir "${params.outdir}/gtdbtk", mode: 'copy', overwrite: true
    
    input:
    path fasta_files
    val  gtdbtk_db
    
    output:
    path "gtdbtk_output/gtdbtk.*.summary.tsv", emit: summary
    path "gtdbtk_output/align",                emit: alignments
    path "gtdbtk_output/*.log",                emit: logs
    
    script:
    """
    set -euo pipefail
    mkdir -p genome_input gtdbtk_output gtdbtk_tmp
    
    for fasta in ${fasta_files}; do
        cp "\$fasta" genome_input/
    done
    
    export GTDBTK_DATA_PATH="${gtdbtk_db}"
    export OMP_NUM_THREADS=${task.cpus}
    
    gtdbtk classify_wf \\
      --genome_dir genome_input \\
      --out_dir gtdbtk_output \\
      --tmpdir gtdbtk_tmp \\
      --cpus ${task.cpus} \\
      --skip_ani_screen \\
      --extension fasta \\
      --debug
    """
}
