#!/usr/bin/env python3
"""
Build the canonical contig annotation master table for the NCBI_4395 batch.

This script collects the scattered table-building steps into one reproducible
entry point. The early table joins are rebuilt from the raw project inputs.
The final CRBC/GTDB-corrected table is treated as a curated handoff because the
repository currently does not contain a deterministic script for that last
taxonomy/type correction step.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Iterable

import pandas as pd


DEFAULT_BATCH = "NCBI_4395_Batch"
LEGACY_RESULT_FINAL_NAME = "Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv"
CANONICAL_FINAL_NAME = "05_master_contig_annotation_table.csv"
MASTER_TABLE_DIR_NAME = "Master_Table"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build and validate the canonical contig annotation master table."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Repository root. Defaults to the directory containing this script.",
    )
    parser.add_argument(
        "--batch",
        default=DEFAULT_BATCH,
        help=f"Batch name under Collect/ and Result/. Default: {DEFAULT_BATCH}",
    )
    parser.add_argument(
        "--data-dir",
        "--output-dir",
        dest="data_dir",
        type=Path,
        default=None,
        help="Directory for normalized inputs, intermediates, final output, and manifest. Defaults to Collect/<batch>/Master_Table.",
    )
    parser.add_argument(
        "--curated-final",
        type=Path,
        default=None,
        help="Curated final CSV used as the authoritative final handoff. Defaults to Collect/<batch>/Master_Table/final/05_master_contig_annotation_table.csv.",
    )
    parser.add_argument(
        "--final-reference",
        type=Path,
        default=None,
        help="Existing final CSV to compare against by MD5.",
    )
    return parser.parse_args()


def repo_path(repo_root: Path, maybe_relative: Path) -> Path:
    return maybe_relative if maybe_relative.is_absolute() else repo_root / maybe_relative


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def line_count(path: Path) -> int:
    with path.open("rb") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def csv_columns(path: Path) -> list[str]:
    return list(pd.read_csv(path, nrows=0).columns)


def file_record(label: str, path: Path) -> dict[str, object]:
    return {
        "label": label,
        "path": str(path),
        "rows": line_count(path),
        "columns": len(csv_columns(path)),
        "md5": md5sum(path),
        "bytes": path.stat().st_size,
    }


def read_csv_text(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def join_non_empty(values: Iterable[str]) -> str:
    return ",".join(value for value in values if value)


def aggregate_by_contig(input_csv: Path, output_csv: Path) -> pd.DataFrame:
    df = read_csv_text(input_csv)
    if "Contig_ID" not in df.columns:
        raise ValueError(f"{input_csv} does not contain a Contig_ID column.")

    grouped = (
        df.groupby("Contig_ID", sort=False, as_index=False)
        .agg(lambda series: join_non_empty(series.astype(str)))
    )
    write_csv(grouped, output_csv)
    return grouped


def merge_mapping_annotations_metadata(
    mapping_csv: Path,
    defense_by_contig_csv: Path,
    amr_by_contig_csv: Path,
    metadata_csv: Path,
    output_csv: Path,
) -> pd.DataFrame:
    mapping = read_csv_text(mapping_csv)
    defense = read_csv_text(defense_by_contig_csv)
    amr = read_csv_text(amr_by_contig_csv)
    metadata = read_csv_text(metadata_csv)

    merged = mapping.merge(defense, on="Contig_ID", how="left")
    merged = merged.merge(amr, on="Contig_ID", how="left")
    merged = merged.merge(metadata, on="Sample_ID", how="left")
    merged = merged.fillna("")
    write_csv(merged, output_csv)
    return merged


def extract_provirus_regions(contig_id: str, contig_type: str) -> list[tuple[int, int]]:
    if "provirus" not in str(contig_type).lower():
        return []
    return [
        (int(start), int(end))
        for start, end in re.findall(r"provirus_region=(\d+)-(\d+)", str(contig_id))
    ]


def parse_positions(value: object) -> list[int]:
    text = str(value).strip()
    if not text:
        return []
    positions = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            positions.append(int(item))
        except ValueError:
            continue
    return positions


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return not (a_end < b_start or b_end < a_start)


def any_feature_overlaps(
    starts: list[int], ends: list[int], regions: list[tuple[int, int]]
) -> bool:
    if not starts or not ends or not regions:
        return False
    for start, end in zip(starts, ends):
        for region_start, region_end in regions:
            if overlaps(start, end, region_start, region_end):
                return True
    return False


def annotate_row_provirus_overlap(row: pd.Series) -> str:
    regions = extract_provirus_regions(row.get("Contig_ID", ""), row.get("Contig_Type", ""))
    checks = [
        ("Defense_on_Provirus", "Defense_Start", "Defense_End"),
        ("PDC_on_Provirus", "PDC_Start", "PDC_End"),
        ("AntiDS_on_Provirus", "AntiDS_Start", "AntiDS_End"),
        ("AMR_on_Provirus", "AMR_Start", "AMR_End"),
    ]
    labels = []
    for label, start_col, end_col in checks:
        if any_feature_overlaps(
            parse_positions(row.get(start_col, "")),
            parse_positions(row.get(end_col, "")),
            regions,
        ):
            labels.append(label)
    return ",".join(labels)


def add_provirus_overlap(input_csv: Path, output_csv: Path) -> pd.DataFrame:
    df = read_csv_text(input_csv)
    df["Provirus_Overlap"] = df.apply(annotate_row_provirus_overlap, axis=1)
    write_csv(df, output_csv)
    return df


def write_checksums(records: list[dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for record in records:
            handle.write(f"{record['md5']}  {record['path']}\n")


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    batch_dir = repo_root / "Collect" / args.batch
    master_table_dir = batch_dir / MASTER_TABLE_DIR_NAME
    data_dir = (
        repo_path(repo_root, args.data_dir).resolve()
        if args.data_dir
        else master_table_dir
    )
    input_dir = data_dir / "input"
    result_dir = repo_root / "Result" / args.batch
    intermediate_dir = data_dir / "intermediate"
    final_dir = data_dir / "final"
    manifest_dir = data_dir / "manifest"
    canonical_final = data_dir / "final" / CANONICAL_FINAL_NAME
    legacy_result_final = result_dir / LEGACY_RESULT_FINAL_NAME

    curated_final = repo_path(
        repo_root,
        args.curated_final
        or (canonical_final if canonical_final.exists() else legacy_result_final),
    ).resolve()
    final_reference = repo_path(
        repo_root,
        args.final_reference
        or curated_final,
    ).resolve()

    input_files = {
        "contig_sample_map": input_dir / "01_contig_sample_map.raw.csv",
        "defense_raw": input_dir / "02_defense_annotations.raw.csv",
        "amr_raw": input_dir / "03_amr_annotations.raw.csv",
        "sample_metadata": input_dir / "04_sample_metadata.raw.csv",
    }

    for label, path in input_files.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing input for {label}: {path}")
    for label, path in {
        "curated_final": curated_final,
        "final_reference": final_reference,
    }.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing {label}: {path}")

    defense_by_contig = intermediate_dir / "01_defense_annotations.by_contig.csv"
    amr_by_contig = intermediate_dir / "02_amr_annotations.by_contig.csv"
    merged_annotations = intermediate_dir / "03_contig_sample_annotations.merged.csv"
    provirus_annotations = intermediate_dir / "04_contig_sample_annotations.provirus_overlap.csv"
    final_table = final_dir / CANONICAL_FINAL_NAME

    print("Step 1/5: aggregating defense annotations by Contig_ID")
    aggregate_by_contig(input_files["defense_raw"], defense_by_contig)

    print("Step 2/5: aggregating AMR annotations by Contig_ID")
    aggregate_by_contig(input_files["amr_raw"], amr_by_contig)

    print("Step 3/5: merging mapping, annotations, and sample metadata")
    merge_mapping_annotations_metadata(
        input_files["contig_sample_map"],
        defense_by_contig,
        amr_by_contig,
        input_files["sample_metadata"],
        merged_annotations,
    )

    print("Step 4/5: annotating feature overlap with provirus regions")
    add_provirus_overlap(merged_annotations, provirus_annotations)

    print("Step 5/5: copying curated CRBC/GTDB-corrected final table")
    final_dir.mkdir(parents=True, exist_ok=True)
    if curated_final != final_table.resolve():
        shutil.copy2(curated_final, final_table)
    else:
        print("Canonical final is already in the output directory; keeping it in place.")

    records = [
        file_record("raw_contig_sample_map", input_files["contig_sample_map"]),
        file_record("raw_defense_annotations", input_files["defense_raw"]),
        file_record("raw_amr_annotations", input_files["amr_raw"]),
        file_record("raw_sample_metadata", input_files["sample_metadata"]),
        file_record("defense_by_contig", defense_by_contig),
        file_record("amr_by_contig", amr_by_contig),
        file_record("merged_annotations", merged_annotations),
        file_record("provirus_annotations", provirus_annotations),
        file_record("curated_final_source", curated_final),
        file_record("normalized_final_output", final_table),
        file_record("canonical_final_reference", final_reference),
    ]
    write_checksums(records, manifest_dir / "checksums.md5")

    manifest = {
        "repo_root": str(repo_root),
        "batch": args.batch,
        "data_dir": str(data_dir),
        "steps": [
            "Read normalized input files from Collect/<batch>/Master_Table/input.",
            "Aggregate 02_defense_annotations.raw.csv by Contig_ID.",
            "Aggregate 03_amr_annotations.raw.csv by Contig_ID.",
            "Merge 01_contig_sample_map.raw.csv, annotation tables, and 04_sample_metadata.raw.csv.",
            "Add Provirus_Overlap annotations.",
            "Use the curated CRBC/GTDB-corrected final table as the authoritative final handoff.",
            "Validate normalized final output against the existing final table by MD5.",
        ],
        "files": records,
    }
    manifest_dir.mkdir(parents=True, exist_ok=True)
    (manifest_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    output_md5 = md5sum(final_table)
    reference_md5 = md5sum(final_reference)
    print(f"Normalized final: {final_table}")
    print(f"Reference final:  {final_reference}")
    print(f"MD5 normalized final: {output_md5}")
    print(f"MD5 reference final:  {reference_md5}")
    if output_md5 != reference_md5:
        raise SystemExit("MD5 validation failed: normalized final does not match reference final.")
    print("MD5 validation passed.")


if __name__ == "__main__":
    main()
