#!/usr/bin/env python3
"""
Compare transferred vs non-transferred plasmid gene sets and draw Venn plots.

This script:
1) Filters rows where Contig_Type3 == Plasmid.
2) Splits plasmid rows into transferred and non-transferred groups.
3) Extracts unique genes for Defense, AMR, and AntiDefense (comma/semicolon separated).
4) Calculates intersections and set sizes.
5) Saves summary tables and two-set Venn-style plots.

Usage example:
  python 223_compare_transferred_plasmid_gene_venn.py \
    -i Result/NCBI_4395_Batch/07_Network/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected_with_Transferred.csv \
    -o Result/NCBI_4395_Batch/07_Network/transferred_plasmid_gene_overlap
"""

import argparse
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle
from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    """Step 1. Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare gene set overlap between transferred and non-transferred plasmids.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Example:\n"
            "  python 223_compare_transferred_plasmid_gene_venn.py \\\n"
            "    -i Result/NCBI_4395_Batch/07_Network/"
            "Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected_with_Transferred.csv \\\n"
            "    -o Result/NCBI_4395_Batch/07_Network/transferred_plasmid_gene_overlap"
        ),
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input CSV file path.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Output directory for statistics and Venn plots.",
    )
    parser.add_argument(
        "--contig-col",
        default="Contig_Type3",
        help="Contig type column name (default: Contig_Type3).",
    )
    parser.add_argument(
        "--contig-value",
        default="Plasmid",
        help="Contig type value to keep (default: Plasmid).",
    )
    parser.add_argument(
        "--transfer-col",
        default="Transferred_Plasmid",
        help="Transfer flag column name (default: Transferred_Plasmid).",
    )
    parser.add_argument(
        "--transfer-yes",
        default="Yes",
        help="Value indicating transferred plasmid (default: Yes).",
    )
    parser.add_argument(
        "--transfer-no",
        default="No",
        help="Value indicating non-transferred plasmid (default: No).",
    )
    parser.add_argument(
        "--defense-col",
        default="Defense_Subtype",
        help="Defense gene column name (default: Defense_Subtype).",
    )
    parser.add_argument(
        "--amr-col",
        default="AMR_Type",
        help="AMR gene column name (default: AMR_Type).",
    )
    parser.add_argument(
        "--antids-col",
        default="AntiDS_Type",
        help="Anti-defense gene column name (default: AntiDS_Type).",
    )
    return parser.parse_args()


def split_gene_tokens(value: str) -> List[str]:
    """Step 2. Split mixed-delimiter gene strings into cleaned tokens."""
    if value is None or pd.isna(value):
        return []
    parts = re.split(r"[;,]", str(value))
    return [p.strip() for p in parts if p.strip()]


def filter_tokens_by_type(tokens: List[str], gene_type: str) -> List[str]:
    """Step 3. Apply gene-type-specific token filters."""
    if gene_type != "Defense":
        return tokens

    exclude_list = {"VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA"}
    return [t for t in tokens if t not in exclude_list and t.lower() != "other"]


def collect_gene_counts(
    series: pd.Series,
    gene_type: str,
) -> Counter:
    """Step 4. Collect gene occurrence counts from a Series."""
    gene_counter: Counter = Counter()

    for value in series:
        tokens = split_gene_tokens(value)
        tokens = filter_tokens_by_type(tokens, gene_type)
        gene_counter.update(tokens)

    return gene_counter


def draw_two_set_venn(
    transferred_set: Set[str],
    non_transferred_set: Set[str],
    title: str,
    out_path: Path,
) -> None:
    """Step 5. Draw a simple two-set Venn-style plot."""
    only_t = len(transferred_set - non_transferred_set)
    only_n = len(non_transferred_set - transferred_set)
    both = len(transferred_set & non_transferred_set)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect("equal")

    circle_left = Circle((0.42, 0.5), 0.28, alpha=0.35, color="#4C72B0")
    circle_right = Circle((0.58, 0.5), 0.28, alpha=0.35, color="#DD8452")
    ax.add_patch(circle_left)
    ax.add_patch(circle_right)

    ax.text(0.31, 0.5, f"{only_t}", ha="center", va="center", fontsize=14)
    ax.text(0.69, 0.5, f"{only_n}", ha="center", va="center", fontsize=14)
    ax.text(0.50, 0.5, f"{both}", ha="center", va="center", fontsize=14, fontweight="bold")

    ax.text(0.30, 0.82, "Transferred", ha="center", va="center", fontsize=12)
    ax.text(0.70, 0.82, "Non-Transferred", ha="center", va="center", fontsize=12)
    ax.set_title(title, fontsize=14)

    ax.set_xlim(0.05, 0.95)
    ax.set_ylim(0.08, 0.92)
    ax.axis("off")

    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def compute_overlap_stats(
    transferred: Set[str],
    non_transferred: Set[str],
) -> Dict[str, int]:
    """Step 6. Compute overlap statistics."""
    inter = transferred & non_transferred
    return {
        "Transferred_Unique_Count": len(transferred - non_transferred),
        "NonTransferred_Unique_Count": len(non_transferred - transferred),
        "Intersection_Count": len(inter),
        "Transferred_Total_Unique": len(transferred),
        "NonTransferred_Total_Unique": len(non_transferred),
        "Union_Count": len(transferred | non_transferred),
    }


def main() -> None:
    """Step 7. Main workflow."""
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, dtype=str, low_memory=False)

    required_cols = [
        args.contig_col,
        args.transfer_col,
        args.defense_col,
        args.amr_col,
        args.antids_col,
    ]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {', '.join(missing_cols)}")

    # Step 8. Filter plasmids and split by transfer status.
    plasmid_df = df[df[args.contig_col] == args.contig_value].copy()
    transferred_df = plasmid_df[plasmid_df[args.transfer_col] == args.transfer_yes].copy()
    non_transferred_df = plasmid_df[plasmid_df[args.transfer_col] == args.transfer_no].copy()

    analyses: Tuple[Tuple[str, str], ...] = (
        ("Defense", args.defense_col),
        ("AMR", args.amr_col),
        ("AntiDefense", args.antids_col),
    )

    summary_rows = []

    for gene_type, col in tqdm(analyses, desc="Analyzing gene types", unit="type"):
        t_counter = collect_gene_counts(transferred_df[col], gene_type)
        n_counter = collect_gene_counts(non_transferred_df[col], gene_type)

        t_set = set(t_counter.keys())
        n_set = set(n_counter.keys())
        stats = compute_overlap_stats(t_set, n_set)

        summary_rows.append({"Gene_Type": gene_type, **stats})

        inter_list = sorted(t_set & n_set)
        only_t_list = sorted(t_set - n_set)
        only_n_list = sorted(n_set - t_set)

        pd.DataFrame(
            {
                "Gene": inter_list,
                "Transferred_Count": [int(t_counter[g]) for g in inter_list],
                "NonTransferred_Count": [int(n_counter[g]) for g in inter_list],
                "Total_Count": [int(t_counter[g] + n_counter[g]) for g in inter_list],
            }
        ).to_csv(
            output_dir / f"{gene_type}_intersection_genes.csv", index=False
        )
        pd.DataFrame(
            {
                "Gene": only_t_list,
                "Transferred_Count": [int(t_counter[g]) for g in only_t_list],
                "Total_Count": [int(t_counter[g]) for g in only_t_list],
            }
        ).to_csv(
            output_dir / f"{gene_type}_transferred_only_genes.csv", index=False
        )
        pd.DataFrame(
            {
                "Gene": only_n_list,
                "NonTransferred_Count": [int(n_counter[g]) for g in only_n_list],
                "Total_Count": [int(n_counter[g]) for g in only_n_list],
            }
        ).to_csv(
            output_dir / f"{gene_type}_non_transferred_only_genes.csv", index=False
        )

        all_genes = sorted(t_set | n_set)
        category_map = {}
        for gene in all_genes:
            if gene in t_set and gene in n_set:
                category_map[gene] = "Intersection"
            elif gene in t_set:
                category_map[gene] = "Transferred_Only"
            else:
                category_map[gene] = "NonTransferred_Only"

        pd.DataFrame(
            {
                "Gene": all_genes,
                "Category": [category_map[g] for g in all_genes],
                "Transferred_Count": [int(t_counter[g]) for g in all_genes],
                "NonTransferred_Count": [int(n_counter[g]) for g in all_genes],
                "Total_Count": [int(t_counter[g] + n_counter[g]) for g in all_genes],
            }
        ).to_csv(output_dir / f"{gene_type}_gene_frequency_table.csv", index=False)

        draw_two_set_venn(
            transferred_set=t_set,
            non_transferred_set=n_set,
            title=f"{gene_type} Gene Overlap on Plasmids",
            out_path=output_dir / f"{gene_type}_venn.png",
        )

    # Step 9. Save summary tables.
    pd.DataFrame(summary_rows).to_csv(
        output_dir / "Transferred_vs_NonTransferred_Plasmid_Gene_Overlap_Summary.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "Metric": [
                "Total_Plasmid_Rows",
                "Transferred_Plasmid_Rows",
                "NonTransferred_Plasmid_Rows",
            ],
            "Value": [
                len(plasmid_df),
                len(transferred_df),
                len(non_transferred_df),
            ],
        }
    ).to_csv(output_dir / "Plasmid_Group_Row_Counts.csv", index=False)

    print(f"Input file: {args.input}")
    print(f"Output directory: {output_dir}")
    print(f"Plasmid rows: {len(plasmid_df)}")
    print(f"Transferred plasmid rows: {len(transferred_df)}")
    print(f"Non-transferred plasmid rows: {len(non_transferred_df)}")


if __name__ == "__main__":
    main()
