process EXTRACT_CONTIG_MAPPING {
    tag "Extract contig mapping from labeled fastas"
    publishDir "${params.outdir}/contig_mapping", mode: 'copy', overwrite: true
    
    input:
    path(labeled_fastas)  // List of all labeled fasta files
    
    output:
    path("contig_mapping.csv"), emit: mapping_table
    
    script:
    """
    # Create temporary directory to hold all labeled fasta files
    mkdir -p labeled_dir
    
    # Link all input files to temporary directory (preserve original filenames)
    for fasta in ${labeled_fastas}; do
        ln -sf "\$(readlink -f \${fasta})" labeled_dir/
    done
    
    # Run Python script to extract mappings
    extract_contig_mapping.py \\
        -i labeled_dir \\
        -o contig_mapping.csv
    """
}

