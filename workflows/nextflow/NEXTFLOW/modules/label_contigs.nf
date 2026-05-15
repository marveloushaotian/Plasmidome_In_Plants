process LABEL_CONTIGS {
    tag "$sample_id"
    publishDir "${params.outdir}/labeled", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(fasta), path(genomad_dir)
    
    output:
    tuple val(sample_id), path("${sample_id}_labeled.fasta"), emit: labeled_fasta
    
    script:
    """
    set -euo pipefail
    
    SUMDIR="${sample_id}_genomad/${sample_id}_filtered_summary/"
    PLASMID_TSV=""
    VIRUS_TSV=""
    
    if [ -d "\$SUMDIR" ]; then
        if compgen -G "\$SUMDIR/*_plasmid_summary.tsv" > /dev/null; then
            for f in "\$SUMDIR"/*_plasmid_summary.tsv; do
              [ -e "\$f" ] && PLASMID_TSV="\$f" && break
            done
        fi
        if compgen -G "\$SUMDIR/*_virus_summary.tsv" > /dev/null; then
            for f in "\$SUMDIR"/*_virus_summary.tsv; do
              [ -e "\$f" ] && VIRUS_TSV="\$f" && break
            done
        fi
    fi
    
    : > plasmid_map.txt
    : > virus_map.txt
    
    if [ -n "\$PLASMID_TSV" ]; then
        awk 'NR>1 {print \$1, \$2}' "\$PLASMID_TSV" > plasmid_map.txt
    fi
    if [ -n "\$VIRUS_TSV" ]; then
        awk 'NR>1 {print \$1, \$2}' "\$VIRUS_TSV" > virus_map.txt
    fi
    
    awk -v pfile="plasmid_map.txt" -v vfile="virus_map.txt" '
    BEGIN {
        while ((getline < pfile) > 0) { plasmid[\$1] = \$2 }
        close(pfile)
        while ((getline < vfile) > 0) {
            split(\$1, a, "|")
            vhit[a[1]] = vhit[a[1]] " " \$1
            vscore[\$1] = \$2
        }
        close(vfile)
    }
    /^>/ {
        cname = substr(\$0, 2)
        label = "chromosome"
        info  = ""
        
        has_p    = (cname in plasmid)
        has_vany = (cname in vhit)
        has_vpro = (has_vany && vhit[cname] ~ /provirus_/)
        
        if (has_p) {
            label = "plasmid"
            info  = "confidence=" plasmid[cname]
        }
        if (has_vany) {
            if (label == "plasmid") { label = "plasmid_virus" } else { label = "virus" }
            if (has_vpro) {
                n = split(vhit[cname], arr, " ")
                for (i=1; i<=n; i++) {
                    if (arr[i] ~ /provirus_/) {
                        split(arr[i], pv, "\\\\|provirus_")
                        split(pv[2], coords, "_")
                        info = info ((info!="" ? "|" : "")) "provirus_region=" coords[1] "-" coords[2]
                    }
                }
            }
        }
        
        print ">" cname "|" label (info != "" ? "|" info : "")
        next
    }
    { print }
    ' "${fasta}" > "${sample_id}_labeled.fasta"
    """
}
