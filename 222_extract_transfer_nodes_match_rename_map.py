#!/usr/bin/env python3
"""
Merge transfer-network source/target IDs, match with rename-map new IDs,
and build a deduplicated original-after-first-pipe list.

Workflow:
1) Read multiple edge TSV files and collect source + target into one column.
2) Save merged list, then save a deduplicated unique list.
3) Match the unique list against rename_map.tsv column "new".
4) Save matched rows from rename_map.tsv.
5) From matched rows, extract "original" after the first "|" and deduplicate.
6) Save the final unique extracted list.

Usage example:
  python 224_extract_transfer_nodes_match_rename_map.py \
    -e Result/NCBI_4395_Batch/07_Network/transfer_network/Medicago_linkage_contig_edges.tsv \
       Result/NCBI_4395_Batch/07_Network/transfer_network/Oryza_linkage_contig_edges.tsv \
       Result/NCBI_4395_Batch/07_Network/transfer_network/Triticum_linkage_contig_edges.tsv \
       Result/NCBI_4395_Batch/07_Network/transfer_network/Zea_linkage_contig_edges.tsv \
    -r Result/NCBI_4395_Batch/07_Network/transfer_network/mmseq_overall_contigs_cluster.rename_map.tsv \
    -o Result/NCBI_4395_Batch/07_Network/transfer_network
"""

import argparse
from pathlib import Path
from typing import Iterable, List, Optional

import pandas as pd
from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    """Step 1. Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge source/target from edge files, match to rename-map 'new', "
            "and export deduplicated original-after-first-pipe list."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-e",
        "--edges",
        nargs="+",
        required=True,
        help="Input edge TSV files (must contain 'source' and 'target').",
    )
    parser.add_argument(
        "-r",
        "--rename-map",
        required=True,
        help="Rename-map TSV file (must contain 'original' and 'new').",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Output directory path.",
    )
    parser.add_argument(
        "--merged-name",
        default="transfer_edges_source_target_merged.tsv",
        help="Output filename for merged source+target list (with duplicates).",
    )
    parser.add_argument(
        "--unique-name",
        default="transfer_edges_source_target_unique.tsv",
        help="Output filename for deduplicated source+target list.",
    )
    parser.add_argument(
        "--matched-name",
        default="mmseq_overall_contigs_cluster.rename_map.matched_by_transfer_edges.tsv",
        help="Output filename for matched rename_map rows.",
    )
    parser.add_argument(
        "--final-list-name",
        default="mmseq_overall_contigs_cluster.original_after_first_pipe.unique.from_matched.tsv",
        help="Output filename for final deduplicated extracted list.",
    )
    return parser.parse_args()


def extract_after_first_pipe(value: str) -> Optional[str]:
    """Step 2. Extract substring after first '|'."""
    if value is None or pd.isna(value):
        return None
    text = str(value)
    idx = text.find("|")
    if idx == -1 or idx == len(text) - 1:
        return None
    return text[idx + 1 :]


def stable_unique(values: Iterable[str]) -> List[str]:
    """Step 3. Deduplicate while preserving first appearance order."""
    seen = set()
    output: List[str] = []
    for v in values:
        if v not in seen:
            seen.add(v)
            output.append(v)
    return output


def main() -> None:
    """Step 4. Execute workflow."""
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    merged_values: List[str] = []

    # Step 5. Read all edge files and collect source + target.
    for edge_file in tqdm(args.edges, desc="Reading edge files", unit="file"):
        edge_df = pd.read_csv(edge_file, sep="\t", dtype=str)
        required = {"source", "target"}
        missing = required - set(edge_df.columns)
        if missing:
            raise ValueError(f"Missing columns in {edge_file}: {', '.join(sorted(missing))}")

        merged_values.extend(edge_df["source"].dropna().astype(str).tolist())
        merged_values.extend(edge_df["target"].dropna().astype(str).tolist())

    # Step 6. Save merged (with duplicates).
    merged_path = output_dir / args.merged_name
    pd.DataFrame({"contig_id": merged_values}).to_csv(merged_path, sep="\t", index=False)

    # Step 7. Deduplicate and save unique list.
    unique_values = stable_unique(merged_values)
    unique_path = output_dir / args.unique_name
    pd.DataFrame({"contig_id": unique_values}).to_csv(unique_path, sep="\t", index=False)

    # Step 8. Match unique IDs with rename_map new column and save matched rows.
    rename_df = pd.read_csv(args.rename_map, sep="\t", dtype=str)
    required_map_cols = {"original", "new"}
    missing_map_cols = required_map_cols - set(rename_df.columns)
    if missing_map_cols:
        raise ValueError(
            f"Missing columns in rename map: {', '.join(sorted(missing_map_cols))}"
        )

    unique_set = set(unique_values)
    matched_df = rename_df[rename_df["new"].astype(str).isin(unique_set)].copy()
    matched_path = output_dir / args.matched_name
    matched_df.to_csv(matched_path, sep="\t", index=False)

    # Step 9. Extract original after first pipe from matched rows and deduplicate.
    extracted_raw = []
    for value in tqdm(matched_df["original"], desc="Extracting original-after-first-pipe", unit="row"):
        extracted = extract_after_first_pipe(value)
        if extracted:
            extracted_raw.append(extracted)

    extracted_unique = stable_unique(extracted_raw)
    final_list_path = output_dir / args.final_list_name
    pd.DataFrame({"extracted_after_first_pipe": extracted_unique}).to_csv(
        final_list_path, sep="\t", index=False
    )

    # Step 10. Print concise summary.
    print(f"Merged rows (source+target, with duplicates): {len(merged_values)}")
    print(f"Unique source+target IDs: {len(unique_values)}")
    print(f"Matched rename_map rows by 'new': {len(matched_df)}")
    print(f"Final unique extracted list size: {len(extracted_unique)}")
    print(f"Saved: {merged_path}")
    print(f"Saved: {unique_path}")
    print(f"Saved: {matched_path}")
    print(f"Saved: {final_list_path}")


if __name__ == "__main__":
    main()
