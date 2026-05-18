#!/usr/bin/env python3
"""
Script name: 101_build_master_table.py

Description:
    Build the canonical contig annotation master table for the NCBI_4395 batch.
    The script rebuilds all table joins from normalized project inputs, adds the
    header-derived GeNomad_Contig_Type, adds the locus-mapping-derived Locus_Mapped_Contig_Type,
    and writes the final master table with the current contig type schema.

Inputs:
    --repo-root         Repository root.
    --batch             Batch name under Collect/ and Result/.
    --data-dir          Master table working directory.
    --locus-contig-map  Unique locus-to-contig type reference table.
    --mobtyper-table    MOB-typer contig mobility table.
    --final-reference   Optional existing final CSV used for validation.

Outputs:
    Master_Table/intermediate/01_defense_annotations.by_contig.csv
    Master_Table/intermediate/02_amr_annotations.by_contig.csv
    Master_Table/intermediate/03_contig_sample_annotations.merged.csv
    Master_Table/intermediate/04_contig_sample_annotations.provirus_overlap.csv
    Master_Table/final/07_contig_annotation_master_table.csv
    Master_Table/manifest/checksums.md5
    Master_Table/manifest/run_manifest.json

Author: Haotian
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


DEFAULT_BATCH = "NCBI_4395_Batch"
CANONICAL_FINAL_NAME = "07_contig_annotation_master_table.csv"
MASTER_TABLE_DIR_NAME = "Master_Table"
TAXONOMY_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
ANNOTATION_NUMERIC_COLUMNS = [
    "Defense_Start",
    "Defense_End",
    "Defense_System_Length",
    "Defense_Gene_Length",
    "Defense_Gene_Num",
    "PDC_Start",
    "PDC_End",
    "PDC_System_Length",
    "PDC_Gene_Length",
    "PDC_Gene_Num",
    "AntiDS_Start",
    "AntiDS_End",
    "AntiDS_Sys_Length",
    "AntiDS_Gene_Length",
    "AntiDS_Gene_Num",
    "AMR_Start",
    "AMR_End",
]


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
        "--locus-contig-map",
        type=Path,
        default=None,
        help="Unique locus-to-contig type reference table. Defaults to Collect/<batch>/Master_Table/input/06_locus_contig_type_reference.raw.tsv.",
    )
    parser.add_argument(
        "--mobtyper-table",
        type=Path,
        default=None,
        help="MOB-typer contig mobility table. Defaults to Collect/<batch>/Master_Table/input/05_mobtyper_contig_mobility.raw.csv.",
    )
    parser.add_argument(
        "--final-reference",
        type=Path,
        default=None,
        help="Existing final CSV to compare against.",
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
    separator = "\t" if path.suffix == ".tsv" else ","
    return list(pd.read_csv(path, sep=separator, nrows=0).columns)


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


def join_annotation_values(series: pd.Series) -> str:
    values = []
    for value in series:
        text = "" if pd.isna(value) else str(value)
        if series.name == "AMR_HMM_Description" and text == "":
            text = "nan"
        if text:
            values.append(text)

    joined = ",".join(values)
    if series.name == "AMR_HMM_Description" and joined == "nan":
        return ""
    return joined


def normalize_legacy_numeric_list(value: str) -> str:
    text = str(value).strip()
    if not text:
        return ""

    parts = [part.strip() for part in text.split(",")]
    if len(parts) == 1:
        return re.sub(r"^(\d+)\.0$", r"\1", parts[0])

    normalized = []
    for part in parts:
        match = re.fullmatch(r"(\d+)(?:\.0)?", part)
        normalized.append(f"{match.group(1)}.0" if match else part)
    return ",".join(normalized)


def normalize_annotation_numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    for column in ANNOTATION_NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = df[column].map(normalize_legacy_numeric_list)
    return df


def normalize_small_float_string(value: str) -> str:
    text = str(value).strip()
    if not text:
        return ""
    try:
        numeric_value = float(text)
    except ValueError:
        return text
    if 0 < abs(numeric_value) < 0.0001:
        return str(numeric_value)
    return text


def aggregate_by_contig(input_csv: Path, output_csv: Path) -> pd.DataFrame:
    df = read_csv_text(input_csv)
    if "Contig_ID" not in df.columns:
        raise ValueError(f"{input_csv} does not contain a Contig_ID column.")

    grouped = (
        df.groupby("Contig_ID", sort=False, as_index=False)
        .agg(join_annotation_values)
    )
    grouped = normalize_annotation_numeric_columns(grouped)
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
    if "virus_percent" in metadata.columns:
        metadata["virus_percent"] = metadata["virus_percent"].map(normalize_small_float_string)

    merged = mapping.merge(defense, on="Contig_ID", how="left")
    merged = merged.merge(amr, on="Contig_ID", how="left")
    merged = merged.merge(metadata, on="Sample_ID", how="left")
    merged = merged.fillna("")
    write_csv(merged, output_csv)
    return merged


def parse_contig_type_from_header(contig_id: str) -> str:
    text = str(contig_id)
    if "|chromosome" in text:
        base_type = "Chromosome"
    elif "|plasmid" in text:
        base_type = "Plasmid"
    elif "|virus" in text:
        base_type = "Virus"
    else:
        base_type = ""

    if "provirus_region=" not in text or not base_type:
        return base_type
    if base_type == "Plasmid":
        return "Plasmid|provirus"
    return f"{base_type}|Provirus"


def add_header_contig_type(df: pd.DataFrame) -> pd.DataFrame:
    if "GeNomad_Contig_Type" in df.columns:
        df = df.drop(columns=["GeNomad_Contig_Type"])
    insert_at = df.columns.get_loc("Contig_ID") + 1
    df.insert(insert_at, "GeNomad_Contig_Type", df["Contig_ID"].map(parse_contig_type_from_header))
    return df


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
    regions = extract_provirus_regions(row.get("Contig_ID", ""), row.get("GeNomad_Contig_Type", ""))
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
    if "GeNomad_Contig_Type" not in df.columns:
        df = add_header_contig_type(df)
    df["Provirus_Overlap"] = df.apply(annotate_row_provirus_overlap, axis=1)
    write_csv(df, output_csv)
    return df


def clean_crbc_genus(genus: str) -> str:
    text = str(genus).strip()
    if not text:
        return "unknown"
    return re.sub(r"_[A-Z]+$", "", text)


def add_crbc_taxonomy(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {f"{level} (gtdb)": f"{level}_CRBC" for level in TAXONOMY_LEVELS}
    missing = [column for column in rename_map if column not in df.columns]
    if missing:
        raise ValueError(f"Missing GTDB taxonomy columns: {missing}")

    df = df.rename(columns=rename_map)
    raw_taxonomy_columns = [level for level in TAXONOMY_LEVELS if level in df.columns]
    df = df.drop(columns=raw_taxonomy_columns)

    insert_at = df.columns.get_loc("Species_CRBC") + 1
    df.insert(insert_at, "Genus_CRBC_Updated", df["Genus_CRBC"].map(clean_crbc_genus))
    return df


def add_predicted_mobility(df: pd.DataFrame, mobtyper_csv: Path) -> pd.DataFrame:
    mobtyper = pd.read_csv(
        mobtyper_csv,
        usecols=["sample_id", "predicted_mobility"],
        dtype=str,
        keep_default_na=False,
    )
    df = df.merge(
        mobtyper,
        left_on="Contig_ID",
        right_on="sample_id",
        how="left",
        sort=False,
    ).drop(columns=["sample_id"])
    df["predicted_mobility"] = df["predicted_mobility"].fillna("")

    if "Provirus_Overlap" in df.columns:
        columns = list(df.columns)
        columns.remove("predicted_mobility")
        columns.insert(columns.index("Provirus_Overlap"), "predicted_mobility")
        df = df[columns]
    return df


def add_locus_mapped_contig_type(
    df: pd.DataFrame,
    locus_contig_map_tsv: Path,
) -> pd.DataFrame:
    mapping = pd.read_csv(
        locus_contig_map_tsv,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    required_columns = {"Downstream_ID", "category"}
    missing = required_columns - set(mapping.columns)
    if missing:
        raise ValueError(f"Missing columns in {locus_contig_map_tsv}: {sorted(missing)}")

    mapping = mapping.rename(columns={"Downstream_ID": "Contig_ID"})
    before_rows = len(df)
    df = df.merge(mapping[["Contig_ID", "category"]], on="Contig_ID", how="inner", sort=False)
    if df.empty:
        raise ValueError("No rows matched the unique locus-to-contig mapping table.")

    df["Locus_Mapped_Contig_Type"] = df["category"].str.capitalize()
    df = df.drop(columns=["category"])
    print(f"  Retained {len(df)} of {before_rows} rows after applying Locus_Mapped_Contig_Type mapping")
    return df


def write_final_master_table(
    provirus_annotations_csv: Path,
    locus_contig_map_tsv: Path,
    mobtyper_csv: Path,
    output_csv: Path,
) -> pd.DataFrame:
    df = read_csv_text(provirus_annotations_csv)
    df = add_header_contig_type(df)
    df = add_crbc_taxonomy(df)
    df = add_predicted_mobility(df, mobtyper_csv)
    df = add_locus_mapped_contig_type(df, locus_contig_map_tsv)
    write_csv(df, output_csv)
    return df


def validate_final_table(
    final_csv: Path,
    reference_csv: Path,
    reference_df: pd.DataFrame | None = None,
) -> None:
    output_df = read_csv_text(final_csv)
    if reference_df is None:
        reference_df = read_csv_text(reference_csv)

    if list(output_df.columns) != list(reference_df.columns):
        raise SystemExit(
            "Final validation failed: output columns do not match the reference."
        )
    if not output_df.equals(reference_df):
        raise SystemExit(
            "Final validation failed: output content does not match the reference."
        )
    print("Final table validation passed.")


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
    intermediate_dir = data_dir / "intermediate"
    final_dir = data_dir / "final"
    manifest_dir = data_dir / "manifest"
    canonical_final = data_dir / "final" / CANONICAL_FINAL_NAME
    locus_contig_map = repo_path(
        repo_root,
        args.locus_contig_map
        or (input_dir / "06_locus_contig_type_reference.raw.tsv"),
    ).resolve()
    mobtyper_table = repo_path(
        repo_root,
        args.mobtyper_table
        or (input_dir / "05_mobtyper_contig_mobility.raw.csv"),
    ).resolve()
    final_reference = repo_path(
        repo_root,
        args.final_reference
        or canonical_final,
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
        "locus_contig_map": locus_contig_map,
        "mobtyper_table": mobtyper_table,
        "final_reference": final_reference,
    }.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing {label}: {path}")

    defense_by_contig = intermediate_dir / "01_defense_annotations.by_contig.csv"
    amr_by_contig = intermediate_dir / "02_amr_annotations.by_contig.csv"
    merged_annotations = intermediate_dir / "03_contig_sample_annotations.merged.csv"
    provirus_annotations = intermediate_dir / "04_contig_sample_annotations.provirus_overlap.csv"
    final_table = final_dir / CANONICAL_FINAL_NAME
    reference_df = read_csv_text(final_reference)

    print("Step 1/6: aggregating defense annotations by Contig_ID")
    aggregate_by_contig(input_files["defense_raw"], defense_by_contig)

    print("Step 2/6: aggregating AMR annotations by Contig_ID")
    aggregate_by_contig(input_files["amr_raw"], amr_by_contig)

    print("Step 3/6: merging mapping, annotations, and sample metadata")
    merge_mapping_annotations_metadata(
        input_files["contig_sample_map"],
        defense_by_contig,
        amr_by_contig,
        input_files["sample_metadata"],
        merged_annotations,
    )

    print("Step 4/6: annotating feature overlap with provirus regions")
    add_provirus_overlap(merged_annotations, provirus_annotations)

    print("Step 5/6: building final master table")
    write_final_master_table(
        provirus_annotations,
        locus_contig_map,
        mobtyper_table,
        final_table,
    )

    print("Step 6/6: validating final table against the reference")
    validate_final_table(final_table, final_reference, reference_df)

    records = [
        file_record("raw_contig_sample_map", input_files["contig_sample_map"]),
        file_record("raw_defense_annotations", input_files["defense_raw"]),
        file_record("raw_amr_annotations", input_files["amr_raw"]),
        file_record("raw_sample_metadata", input_files["sample_metadata"]),
        file_record("defense_by_contig", defense_by_contig),
        file_record("amr_by_contig", amr_by_contig),
        file_record("merged_annotations", merged_annotations),
        file_record("provirus_annotations", provirus_annotations),
        file_record("unique_locus_contig_map", locus_contig_map),
        file_record("merged_mobtyper_table", mobtyper_table),
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
            "Parse GeNomad_Contig_Type from the labeled Contig_ID header.",
            "Add Provirus_Overlap annotations.",
            "Rename GTDB taxonomy columns to CRBC taxonomy columns and clean Genus_CRBC_Updated.",
            "Join predicted_mobility from input/05_mobtyper_contig_mobility.raw.csv.",
            "Join category from input/06_locus_contig_type_reference.raw.tsv as Locus_Mapped_Contig_Type.",
            "Write the current final schema with GeNomad_Contig_Type and Locus_Mapped_Contig_Type.",
            "Validate normalized final output against the existing final table.",
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
    if output_md5 == reference_md5:
        print("MD5 validation passed.")
    else:
        print("MD5 differs: normalized final does not byte-match the reference.")


if __name__ == "__main__":
    main()
