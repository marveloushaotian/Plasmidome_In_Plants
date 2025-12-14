process IQTREE_GROUPED {
    tag "iqtree_${group_name}_${alignment.simpleName}"
    publishDir "${params.outdir}/trees_grouped/${group_name}", mode: 'copy', overwrite: true

    input:
    tuple val(group_name), path(alignment)

    output:
    path "${alignment.simpleName}.treefile", emit: tree
    path "${alignment.simpleName}.iqtree",   emit: report
    path "${alignment.simpleName}.log",      emit: log

    script:
    """
    set -euo pipefail

    NSEQ=\$(grep -c '^>' "${alignment}" || echo 0)
    if [ "\$NSEQ" -lt 2 ]; then
      echo "[ERROR] ${alignment} has only \$NSEQ sequences; need >=2 for a meaningful tree." >&2
      exit 3
    fi

    echo "[INFO] Building tree for group: ${group_name}"
    echo "[INFO] Alignment file: ${alignment}"
    echo "[INFO] Number of sequences: \$NSEQ"

    if [ ! -s best_model.txt ]; then
      iqtree2 \\
        -s "${alignment}" \\
        -st AA \\
        -m MFP \\
        -pre model_test \\
        -T ${task.cpus} \\
        -mem ${task.memory.toMega()}M
      awk -F':' '/Best-fit model/ {gsub(/[[:space:]]+/,"",\$2); print \$2}' model_test.iqtree > best_model.txt || true
    fi

    BEST_MODEL="LG+R8"
    if [ -s best_model.txt ]; then
      BEST_MODEL=\$(cat best_model.txt)
    fi
    echo "[INFO] Using model: \$BEST_MODEL"

    iqtree2 \\
      -s "${alignment}" \\
      -st AA \\
      -m "\$BEST_MODEL" \\
      -B ${params.bootstrap_reps} \\
      -alrt ${params.alrt_reps} \\
      -bnni \\
      -T ${task.cpus} \\
      -mem ${task.memory.toMega()}M \\
      -pre "${alignment.simpleName}"
    """
}
