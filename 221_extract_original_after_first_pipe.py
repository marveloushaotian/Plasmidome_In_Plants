#!/usr/bin/env python3
"""
Extract unique substrings after the first pipe from a target TSV column.

This script reads a TSV file, takes values from a target column (default:
"original"), removes the first "|" and everything before it, de-duplicates
the extracted strings while preserving first appearance order, and writes
them to a new TSV file.

Usage example:
  python 222_extract_original_after_first_pipe.py \
    -i Result/NCBI_4395_Batch/07_Network/transfer_network/mmseq_overall_contigs_cluster.rename_map.tsv \
    -o Result/NCBI_4395_Batch/07_Network/transfer_network/mmseq_overall_contigs_cluster.original_after_first_pipe.unique.tsv
"""

import argparse
import csv
from pathlib import Path
from typing import Optional

from tqdm import tqdm


def extract_after_first_pipe(text: str) -> Optional[str]:
    """Return substring after the first '|' or None when no pipe exists."""
    if text is None:
        return None
    idx = text.find("|")
    if idx == -1 or idx == len(text) - 1:
        return None
    return text[idx + 1 :]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Extract unique strings after the first '|' from a target "
            "column in a TSV file."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=(
            "Example:\n"
            "  python 222_extract_original_after_first_pipe.py \\\n"
            "    -i Result/NCBI_4395_Batch/07_Network/transfer_network/"
            "mmseq_overall_contigs_cluster.rename_map.tsv \\\n"
            "    -o Result/NCBI_4395_Batch/07_Network/transfer_network/"
            "mmseq_overall_contigs_cluster.original_after_first_pipe.unique.tsv"
        ),
    )
    parser.add_argument("-i", "--input", required=True, help="Input TSV file path.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file path.")
    parser.add_argument(
        "-c",
        "--column",
        default="original",
        help="Target column name to process (default: original).",
    )
    return parser.parse_args()


def main() -> None:
    """Run extraction and save unique values."""
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    unique_values = []
    invalid_rows = 0

    with input_path.open("r", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        if reader.fieldnames is None or args.column not in reader.fieldnames:
            raise ValueError(
                f"Column '{args.column}' not found. Available columns: {reader.fieldnames}"
            )

        for row in tqdm(reader, desc="Extracting", unit="row"):
            extracted = extract_after_first_pipe(row.get(args.column, ""))
            if extracted is None or extracted == "":
                invalid_rows += 1
                continue
            if extracted not in seen:
                seen.add(extracted)
                unique_values.append(extracted)

    with output_path.open("w", encoding="utf-8", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["extracted_after_first_pipe"])
        for value in unique_values:
            writer.writerow([value])

    print(f"Input file: {input_path}")
    print(f"Output file: {output_path}")
    print(f"Unique extracted values: {len(unique_values)}")
    print(f"Skipped invalid rows: {invalid_rows}")


if __name__ == "__main__":
    main()
