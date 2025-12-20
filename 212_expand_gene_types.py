#!/usr/bin/env python3
"""
Expand gene type columns into multiple columns with counts.
Supports: AMR_Type, Defense_Subtype, AntiDS_Type

Each unique gene type becomes a separate column with count values.

Usage:
    python 212_expand_gene_types.py -i <input.csv> -o <output.csv> -t <type>

Arguments:
    -i: Input CSV file path
    -o: Output CSV file path
    -t: Gene type to expand: 'amr', 'defense', or 'antidefense'

Example:
    python 212_expand_gene_types.py -i Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv -o Result/NCBI_4395_Batch/06_Cluter/Contig_Sample_Mapping_Expanded_AMR.csv -t amr
    python 212_expand_gene_types.py -i Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv -o Result/NCBI_4395_Batch/06_Cluter/Contig_Sample_Mapping_Expanded_Defense.csv -t defense
    python 212_expand_gene_types.py -i Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv -o Result/NCBI_4395_Batch/06_Cluter/Contig_Sample_Mapping_Expanded_AntiDefense.csv -t antidefense
"""

import argparse
import pandas as pd
from collections import Counter
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Expand gene type columns into multiple columns with counts.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV file path")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file path")
    parser.add_argument(
        "-t", "--type",
        required=True,
        choices=['amr', 'defense', 'antidefense'],
        help="Gene type to expand: 'amr', 'defense', or 'antidefense'"
    )
    return parser.parse_args()


def count_gene_types(gene_str):
    """
    Parse gene type string and count each type.
    E.g., "tetracycline,tetracycline" -> {"tetracycline": 2}
          "beta-lactam,aminoglycoside" -> {"beta-lactam": 1, "aminoglycoside": 1}
    """
    if pd.isna(gene_str) or gene_str == "":
        return Counter()
    # Split by comma and strip whitespace
    gene_types = [s.strip() for s in str(gene_str).split(",")]
    # Filter out empty strings
    gene_types = [s for s in gene_types if s]
    return Counter(gene_types)


def main():
    args = parse_args()
    
    # Map type to column name and prefix
    type_config = {
        'amr': {
            'column': 'AMR_Type',
            'prefix': 'AMR_',
            'name': 'AMR'
        },
        'defense': {
            'column': 'Defense_Subtype',
            'prefix': 'Defense_',
            'name': 'Defense'
        },
        'antidefense': {
            'column': 'AntiDS_Type',
            'prefix': 'AntiDS_',
            'name': 'Anti-Defense'
        }
    }
    
    config = type_config[args.type]
    column_name = config['column']
    prefix = config['prefix']
    type_name = config['name']
    
    # 1. Read input CSV
    print(f"Reading input file: {args.input}")
    df = pd.read_csv(args.input, low_memory=False)
    print(f"Total rows: {len(df)}")
    
    # Check if column exists
    if column_name not in df.columns:
        print(f"Error: Column '{column_name}' not found in input file")
        return
    
    # 2. Get all unique gene types
    print(f"Extracting all unique {type_name} types...")
    all_gene_types = set()
    for gene_str in tqdm(df[column_name], desc=f"Scanning {type_name} types"):
        counts = count_gene_types(gene_str)
        all_gene_types.update(counts.keys())
    
    all_gene_types = sorted(all_gene_types)
    print(f"Found {len(all_gene_types)} unique {type_name} types")
    
    # 3. Create new columns for each gene type
    print(f"Expanding {column_name} into columns...")
    
    # Initialize new columns with zeros
    for gene_type in all_gene_types:
        df[f"{prefix}{gene_type}"] = 0
    
    # Fill in counts for each row
    for idx in tqdm(range(len(df)), desc="Processing rows"):
        gene_str = df.loc[idx, column_name]
        counts = count_gene_types(gene_str)
        for gene_type, count in counts.items():
            df.loc[idx, f"{prefix}{gene_type}"] = count
    
    # 4. Save output
    print(f"Saving output to: {args.output}")
    df.to_csv(args.output, index=False)
    print(f"Done! Output has {len(df)} rows and {len(df.columns)} columns")
    print(f"Added {len(all_gene_types)} new {type_name} type columns")


if __name__ == "__main__":
    main()

