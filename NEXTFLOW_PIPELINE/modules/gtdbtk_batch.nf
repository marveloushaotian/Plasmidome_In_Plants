// GTDB-Tk batch processing module
// Splits genomes into batches and processes them in parallel to work around pplacer single-thread bottleneck

process GTDBTK_SPLIT {
    tag "split-genomes"
    
    cpus 2
    memory '8 GB'
    conda params.envs.base
    
    input:
    path fasta_files
    
    output:
    path "batch_*", emit: batches
    
    script:
    def batch_size = 50
    """
    #!/bin/bash
    set -euo pipefail
    
    # Collect all fasta files in current directory (staged by Nextflow)
    batch_size=${batch_size}
    batch_num=1
    file_count=0
    
    # Create first batch directory
    mkdir -p batch_\$(printf "%03d" \$batch_num)
    
    # Iterate through all fasta files
    for fasta in *.fasta; do
        if [ -f "\$fasta" ]; then
            # Copy to current batch directory
            cp "\$fasta" batch_\$(printf "%03d" \$batch_num)/
            
            file_count=\$((file_count + 1))
            
            # If batch is full, start new batch
            if [ \$file_count -eq \$batch_size ]; then
                batch_num=\$((batch_num + 1))
                mkdir -p batch_\$(printf "%03d" \$batch_num)
                file_count=0
            fi
        fi
    done
    
    # Remove empty last batch if exists
    if [ \$file_count -eq 0 ] && [ \$batch_num -gt 1 ]; then
        rmdir batch_\$(printf "%03d" \$batch_num) 2>/dev/null || true
    fi
    
    echo "Split genomes into \$batch_num batches (${batch_size} genomes per batch)"
    """
}

process GTDBTK_CLASSIFY_BATCH {
    tag "batch-${batch_dir.name}"
    publishDir "${params.outdir}/gtdbtk_batches/${batch_dir.name}", mode: 'copy', overwrite: true
    
    cpus 4
    memory '40 GB'
    conda params.envs.gtdbtk
    
    input:
    path batch_dir
    val gtdbtk_db
    
    output:
    path "gtdbtk_output/gtdbtk.*.summary.tsv", emit: summary, optional: true
    path "gtdbtk_output/align/*.fasta*", emit: alignments, optional: true
    path "gtdbtk_output/*.log", emit: logs
    
    script:
    """
    set -euo pipefail
    mkdir -p gtdbtk_output gtdbtk_tmp
    
    export GTDBTK_DATA_PATH="${gtdbtk_db}"
    export OMP_NUM_THREADS=${task.cpus}
    
    gtdbtk classify_wf \\
      --genome_dir ${batch_dir} \\
      --out_dir gtdbtk_output \\
      --tmpdir gtdbtk_tmp \\
      --cpus ${task.cpus} \\
      --pplacer_cpus 1 \\
      --skip_ani_screen \\
      --extension fasta \\
      --debug
    """
}

process GTDBTK_MERGE {
    tag "merge-all-batches"
    publishDir "${params.outdir}/gtdbtk", mode: 'copy', overwrite: true
    
    cpus 4
    memory '32 GB'
    conda params.envs.base
    
    input:
    path summary_files
    path alignment_files
    
    output:
    path "gtdbtk.*.summary.tsv", emit: summary
    path "align", emit: alignments
    path "*.log", emit: logs
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import glob
    import os
    import shutil
    
    # Merge bacteria summary files
    bac_files = glob.glob("**/gtdbtk.bac120.summary.tsv", recursive=True)
    if bac_files:
        dfs = [pd.read_csv(f, sep='\\t') for f in bac_files]
        merged = pd.concat(dfs, ignore_index=True)
        merged.to_csv('gtdbtk.bac120.summary.tsv', sep='\\t', index=False)
        print(f"Merged {len(bac_files)} bacteria batches: {len(merged)} genomes total")
    
    # Merge archaea summary files
    ar_files = glob.glob("**/gtdbtk.ar53.summary.tsv", recursive=True)
    if ar_files:
        dfs = [pd.read_csv(f, sep='\\t') for f in ar_files]
        merged = pd.concat(dfs, ignore_index=True)
        merged.to_csv('gtdbtk.ar53.summary.tsv', sep='\\t', index=False)
        print(f"Merged {len(ar_files)} archaea batches: {len(merged)} genomes total")
    
    # Merge alignment files
    os.makedirs('align', exist_ok=True)
    for align_file in glob.glob("**/*.fasta*", recursive=True):
        if '/align/' in align_file or '\\\\align\\\\' in align_file:
            basename = os.path.basename(align_file)
            dest = os.path.join('align', basename)
            
            if dest.endswith('.gz'):
                # For gzipped files, concatenate
                if os.path.exists(dest):
                    os.system(f"cat {align_file} >> {dest}")
                else:
                    shutil.copy(align_file, dest)
            else:
                # For plain fasta, append sequences
                if os.path.exists(dest):
                    with open(dest, 'a') as out_f:
                        with open(align_file, 'r') as in_f:
                            out_f.write(in_f.read())
                else:
                    shutil.copy(align_file, dest)
    
    # Create a summary log
    with open('gtdbtk_batch_merge.log', 'w') as f:
        f.write(f"Bacteria batches merged: {len(bac_files)}\\n")
        f.write(f"Archaea batches merged: {len(ar_files)}\\n")
        f.write(f"Alignment files merged: {len(glob.glob('align/*'))}\\n")
    
    print("Merge completed successfully")
    """
}

