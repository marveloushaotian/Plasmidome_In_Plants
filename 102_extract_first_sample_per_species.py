#!/usr/bin/env python3
"""
Extract first Sample_ID for each unique species

This script reads a CSV file containing contig-sample mapping data and extracts
the first Sample_ID for each unique species value.

Usage:
    python 102_extract_first_sample_per_species.py
"""

import argparse
import pandas as pd
from tqdm import tqdm

DEFAULT_INPUT = "Collect/NCBI_4395_Batch/Master_Table/final/07_contig_annotation_master_table.csv"
DEFAULT_OUTPUT = "Result/NCBI_4395_Batch/05_Tree/Genus_Level/First_Sample_Per_Species.csv"
DEFAULT_SPECIES_COLUMN = "Species_CRBC"


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Extract first Sample_ID for each unique species'
    )
    parser.add_argument(
        '-i', '--input',
        default=DEFAULT_INPUT,
        help='Input CSV file path'
    )
    parser.add_argument(
        '-o', '--output',
        default=DEFAULT_OUTPUT,
        help='Output CSV file path with Sample_ID and species columns'
    )
    parser.add_argument(
        '--species-column',
        default=DEFAULT_SPECIES_COLUMN,
        help='Species column to use (default: Species_CRBC)'
    )
    return parser.parse_args()


def main():
    """Main function"""
    args = parse_args()
    
    # Step 1: Read the CSV file
    print(f"Reading input file: {args.input}")
    df = pd.read_csv(args.input)
    
    print(f"Total rows: {len(df)}")
    
    species_col = args.species_column
    if species_col not in df.columns:
        raise ValueError(f"Column '{species_col}' not found in input file")

    # Step 2: Get unique species values
    print(f"Processing unique species from column: {species_col}")
    
    # Drop rows where species is NaN or empty
    df_filtered = df[df[species_col].notna()]
    
    # Get first Sample_ID for each unique species
    # Using drop_duplicates to keep first occurrence
    unique_species_samples = df_filtered.drop_duplicates(
        subset=[species_col],
        keep='first'
    )[['Sample_ID', species_col]]
    
    print(f"Found {len(unique_species_samples)} unique species")
    
    # Step 3: Save to output CSV file
    print(f"Saving results to: {args.output}")
    unique_species_samples.to_csv(args.output, index=False)
    
    print(f"Done! Saved {len(unique_species_samples)} entries to {args.output}")


if __name__ == '__main__':
    main()
