process EGGNOG_DIAMOND {
    tag "EggNOG diamond search"
    publishDir "${params.outdir}/eggnog/diamond_results", mode: 'copy', overwrite: true
    conda "/scratch/miniforge3/envs/eggnog"
    
    input:
    path(faa_files)  // List of all prokka .faa files
    val(database_dir)
    
    output:
    path("diamond_results"), emit: diamond_dir
    
    script:
    """
    set -e
    
    # 1. Merge all .faa files into one
    echo "Merging ${faa_files.size()} faa files..."
    cat ${faa_files} > merged_proteins.faa
    
    TOTAL_SEQS=\$(grep -c "^>" merged_proteins.faa || true)
    echo "Total sequences: \$TOTAL_SEQS"
    
    # 2. Split using seqkit (activate seqkit environment)
    echo "Splitting into chunks of 200000 sequences..."
    mkdir -p chunks
    
    # Use seqkit from conda environment
    conda run -n seqkit seqkit split2 -s 200000 -O chunks merged_proteins.faa
    
    # Rename chunk files to have .faa extension
    cd chunks
    for f in merged_proteins.part_*.fasta; do
        [ -e "\$f" ] || continue
        mv "\$f" "\${f%.fasta}.faa"
    done
    cd ..
    
    CHUNK_COUNT=\$(ls chunks/*.faa 2>/dev/null | wc -l)
    echo "Created \$CHUNK_COUNT chunk files"
    
    # 3. Run diamond search on all chunks in parallel
    echo "Running diamond search on \$CHUNK_COUNT chunks..."
    mkdir -p diamond_results
    
    # Copy chunks to diamond_results directory
    cp chunks/*.faa diamond_results/
    
    # Run diamond search script
    bash ${projectDir}/bin/run_eggnog_diamond.sh \\
        --chunk_dir diamond_results \\
        --database ${database_dir} \\
        --max_parallel 10 \\
        --cpu 20
    
    echo "Diamond search completed"
    
    # Verify output files
    RESULT_COUNT=\$(ls diamond_results/*.emapper.seed_orthologs 2>/dev/null | wc -l)
    echo "Generated \$RESULT_COUNT diamond result files"
    
    if [ \$RESULT_COUNT -eq 0 ]; then
        echo "Error: No diamond results generated"
        exit 1
    fi
    """
}

process EGGNOG_ANNOTATION {
    tag "EggNOG annotation"
    publishDir "${params.outdir}/eggnog", mode: 'copy', overwrite: true
    conda "/scratch/miniforge3/envs/eggnog"
    
    input:
    path(diamond_dir)
    val(database_dir)
    
    output:
    path("annotation_results"), emit: annotations
    
    script:
    """
    set -e
    
    echo "Starting EggNOG annotation..."
    
    # Count input files
    INPUT_COUNT=\$(ls ${diamond_dir}/*.emapper.seed_orthologs 2>/dev/null | wc -l)
    echo "Found \$INPUT_COUNT diamond result files to annotate"
    
    if [ \$INPUT_COUNT -eq 0 ]; then
        echo "Error: No diamond result files found"
        exit 1
    fi
    
    # Create output directory
    mkdir -p annotation_results
    
    # Run annotation script
    bash ${projectDir}/bin/run_eggnog_annotation.sh \\
        --input_dir ${diamond_dir} \\
        --output_dir annotation_results \\
        --database ${database_dir} \\
        --cpu_per_task 10 \\
        --total_cpu 120
    
    echo "Annotation completed"
    
    # Verify output files
    RESULT_COUNT=\$(ls annotation_results/*_eggnog.emapper.annotations 2>/dev/null | wc -l)
    echo "Generated \$RESULT_COUNT annotation result files"
    
    if [ \$RESULT_COUNT -eq 0 ]; then
        echo "Warning: No annotation results generated, check logs"
    fi
    
    # Summary
    echo "=== Summary ==="
    echo "Input files: \$INPUT_COUNT"
    echo "Output files: \$RESULT_COUNT"
    echo "Log files: \$(ls annotation_results/*.log 2>/dev/null | wc -l)"
    """
}

