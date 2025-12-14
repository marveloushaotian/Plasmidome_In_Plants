process PREP_CHECKM_INPUTS {
    tag "prep-checkm-inputs"
    publishDir "${params.outdir}/checkm_inputs", mode: 'copy', overwrite: true
    
    input:
    path labeled_fastas
    
    output:
    path "checkm_input_all",  emit: bins_all_dir
    path "checkm_input_chro", emit: bins_chro_dir
    
    script:
    """
    set -euo pipefail
    mkdir -p checkm_input_all checkm_input_chro
    
    for f in ${labeled_fastas}; do
      bname=\$(basename "\$f")
      sid=\${bname%_labeled.fasta}
      cp "\$f" "checkm_input_all/\${sid}.fasta"
      
      awk '/^>/{keep=(\$0 ~ /\\\\|chromosome/)} keep' "\$f" > "checkm_input_chro/\${sid}.fasta" || true
      
      if [ ! -s "checkm_input_chro/\${sid}.fasta" ]; then
        rm -f "checkm_input_chro/\${sid}.fasta"
      fi
    done
    
    date -u > "checkm_input_all/_prepared.timestamp"
    date -u > "checkm_input_chro/_prepared.timestamp"
    """
}
