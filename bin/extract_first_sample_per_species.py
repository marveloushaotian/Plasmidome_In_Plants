#!/usr/bin/env python3
"""
Extract first Sample_ID for each unique Species (gtdb)

This script reads a CSV file containing contig-sample mapping data and extracts
the first Sample_ID for each unique Species (gtdb) value.

Usage:
    python extract_first_sample_per_species.py -i input.csv -o output.csv
"""

import argparse
import pandas as pd
from tqdm import tqdm


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Extract first Sample_ID for each unique Species (gtdb)'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CSV file path'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output CSV file path with Sample_ID and Species (gtdb) columns'
    )
    return parser.parse_args()


def main():
    """Main function"""
    args = parse_args()
    
    # Step 1: Read the CSV file
    print(f"Reading input file: {args.input}")
    df = pd.read_csv(args.input)
    
    print(f"Total rows: {len(df)}")
    
    # Step 2: Get unique Species (gtdb) values
    print("Processing unique Species (gtdb)...")
    
    # Drop rows where Species (gtdb) is NaN or empty
    df_filtered = df[df['Species (gtdb)'].notna()]
    
    # Get first Sample_ID for each unique Species (gtdb)
    # Using drop_duplicates to keep first occurrence
    unique_species_samples = df_filtered.drop_duplicates(
        subset=['Species (gtdb)'], 
        keep='first'
    )[['Sample_ID', 'Species (gtdb)']]
    
    print(f"Found {len(unique_species_samples)} unique Species (gtdb)")
    
    # Step 3: Save to output CSV file
    print(f"Saving results to: {args.output}")
    unique_species_samples.to_csv(args.output, index=False)
    
    print(f"Done! Saved {len(unique_species_samples)} entries to {args.output}")


if __name__ == '__main__':
    main()

