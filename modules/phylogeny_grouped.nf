process PHYLOGENY_GROUPED {
    tag "phylogeny_${group_name}"
    publishDir "${params.outdir}/phylogeny_grouped/${group_name}", mode: 'copy', overwrite: true

    input:
    path align_dir
    path hq_list
    tuple val(group_name), val(genome_ids)

    output:
    tuple val(group_name), path("filtered_alignments/*.fasta"), emit: filtered_alignments
    path "alignment_stats.txt",        emit: stats
    path "filter_report.log",          emit: report

    script:
    """
    set -euo pipefail
    mkdir -p filtered_alignments

    ALIGN_DIR='${align_dir}'
    HQ_LIST='${hq_list}'
    GROUP_NAME='${group_name}'
    GENOME_IDS='${genome_ids.join('\n')}'

    echo "[INFO] Processing group: \$GROUP_NAME"
    echo "[INFO] Group members: \$GENOME_IDS"

    pick_first() {
      for f in "\$@"; do
        [ -f "\$f" ] && { printf '%s' "\$f"; return 0; }
      done
      return 1
    }

    BAC_CAND1="\${ALIGN_DIR}/gtdbtk.bac120.user_msa.fasta.gz"
    BAC_CAND2="\${ALIGN_DIR}/gtdbtk.bac120.msa.fasta.gz"
    BAC_CAND3="\${ALIGN_DIR}/gtdbtk.bac120.user_msa.fasta"
    BAC_CAND4="\${ALIGN_DIR}/gtdbtk.bac120.msa.fasta"
    BAC_MSA="\$(pick_first "\$BAC_CAND1" "\$BAC_CAND2" "\$BAC_CAND3" "\$BAC_CAND4" || true)"

    echo "Creating group-specific HQ list for \$GROUP_NAME"
    echo "\$GENOME_IDS" | tr ' ' '\\n' > "hq_list_\${GROUP_NAME}.txt"

    if [ -n "\$BAC_MSA" ]; then
      python3 - "\$BAC_MSA" "hq_list_\${GROUP_NAME}.txt" "filtered_alignments/bacterial_\${GROUP_NAME}.fasta" "filter_report_\${GROUP_NAME}.log" "\$GROUP_NAME" <<'PY'
import sys, gzip, re, os

msa_fp, idlist_fp, out_fp, rpt_fp, group_name = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def header_id(h):
    h = h[1:].strip()
    return re.split(r"\s|\t", h, maxsplit=1)[0]

targets, seen = [], set()
with open(idlist_fp) as f:
    for line in f:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s in seen:
            continue
        seen.add(s)
        targets.append(s)

seqs = {}
with open_any(msa_fp) as f:
    cur_h, buf = None, []
    for line in f:
        if line.startswith(">"):
            if cur_h is not None:
                seqs[cur_h] = "".join(buf)
            cur_h = header_id(line)
            buf = []
        else:
            buf.append(line.strip())
    if cur_h is not None:
        seqs[cur_h] = "".join(buf)

found, missing = 0, []
with open(out_fp, "w") as out:
    for tid in targets:
        s = seqs.get(tid)
        if s is None:
            missing.append(tid)
        else:
            out.write(f">{tid}\\n")
            for i in range(0, len(s), 80):
                out.write(s[i:i+80] + "\\n")
            found += 1

extras = sorted(set(seqs.keys()) - set(targets))
with open(rpt_fp, "w") as log:
    log.write(f"[filter] Group: ${sys.argv[5]}\\n")
    log.write(f"[filter] Targets in list: {len(targets)}\\n")
    log.write(f"[filter] Found in MSA   : {found}\\n")
    log.write(f"[filter] Missing       : {len(missing)}\\n")
    if missing:
        log.write(f"[filter] Missing IDs (first 20): {missing[:20]}\\n")
    log.write(f"[filter] Extra IDs in MSA not in list: {len(extras)}\\n")

if found == 0:
    sys.stderr.write(f"[ERROR] No sequences found for group {group_name}\\n")
    sys.exit(2)
PY
    fi

    ARC_CAND1="\${ALIGN_DIR}/gtdbtk.ar53.user_msa.fasta.gz"
    ARC_CAND2="\${ALIGN_DIR}/gtdbtk.ar53.msa.fasta.gz"
    ARC_CAND3="\${ALIGN_DIR}/gtdbtk.ar53.user_msa.fasta"
    ARC_CAND4="\${ALIGN_DIR}/gtdbtk.ar53.msa.fasta"
    ARC_MSA="\$(pick_first "\$ARC_CAND1" "\$ARC_CAND2" "\$ARC_CAND3" "\$ARC_CAND4" || true)"

    if [ -n "\$ARC_MSA" ]; then
      python3 - "\$ARC_MSA" "hq_list_\${GROUP_NAME}.txt" "filtered_alignments/archaeal_\${GROUP_NAME}.fasta" /dev/null <<'PY'
import sys, gzip, re

msa_fp, idlist_fp, out_fp = sys.argv[1], sys.argv[2], sys.argv[3]
def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")
def header_id(h):
    h = h[1:].strip()
    return re.split(r"\s|\t", h, maxsplit=1)[0]
targets = [s.strip() for s in open(idlist_fp) if s.strip() and not s.startswith("#")]
targets_set = set(targets)
with open_any(msa_fp) as f, open(out_fp, "w") as out:
    keep = False
    for line in f:
        if line.startswith(">"):
            keep = (header_id(line) in targets_set)
            if keep:
                out.write(line)
        else:
            if keep:
                # 80-col wrap
                s = line.strip()
                for i in range(0, len(s), 80):
                    out.write(s[i:i+80] + "\\n")
PY
    fi

    echo "Alignment Statistics for Group: \$GROUP_NAME" > alignment_stats.txt
    found=0
    for f in filtered_alignments/*.fasta; do
      if [ -f "\$f" ]; then
        found=1
        printf '%s: %d sequences\\n' "\$(basename "\$f")" "\$(grep -c '^>' "\$f")" >> alignment_stats.txt
      fi
    done
    if [ "\$found" -eq 0 ]; then
      echo "No filtered alignments were generated for group \$GROUP_NAME." >> alignment_stats.txt
      exit 2
    fi
    """
}
