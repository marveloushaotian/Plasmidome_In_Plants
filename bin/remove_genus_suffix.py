#!/usr/bin/env python3
"""
Remove underscore and everything after it from Genus_CRBC_Updated column

This script processes a CSV file and removes the underscore and any characters
following it in the Genus_CRBC_Updated column (e.g., "Paenibacillus_A" becomes "Paenibacillus").

Usage:
    python remove_genus_suffix.py -i input.csv -o output.csv
"""

import argparse
import pandas as pd
from tqdm import tqdm


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Remove underscore suffix from Genus_CRBC_Updated column'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CSV file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output CSV file'
    )
    return parser.parse_args()


def remove_suffix(genus_name):
    """
    Remove underscore and everything after it from genus name
    
    Args:
        genus_name: String like "Paenibacillus_A" or "Nocardioides"
    
    Returns:
        String with suffix removed (e.g., "Paenibacillus")
    """
    if pd.isna(genus_name):
        return genus_name
    
    genus_str = str(genus_name)
    
    # Split by underscore and take only the first part
    if '_' in genus_str:
        return genus_str.split('_')[0]
    else:
        return genus_str


def main():
    """Main function"""
    args = parse_args()
    
    # Step 1: Read input CSV
    print(f"Reading input file: {args.input}")
    df = pd.read_csv(args.input)
    print(f"Total rows: {len(df)}")
    
    # Check if column exists
    if 'Genus_CRBC_Updated' not in df.columns:
        print("Error: Column 'Genus_CRBC_Updated' not found in the file!")
        print(f"Available columns: {', '.join(df.columns)}")
        return
    
    # Step 2: Count entries with underscore before modification
    original_with_underscore = df['Genus_CRBC_Updated'].astype(str).str.contains('_').sum()
    print(f"Entries with underscore in Genus_CRBC_Updated: {original_with_underscore}")
    
    # Step 3: Process the column
    print("Removing suffixes from Genus_CRBC_Updated column...")
    tqdm.pandas(desc="Processing")
    df['Genus_CRBC_Updated'] = df['Genus_CRBC_Updated'].progress_apply(remove_suffix)
    
    # Step 4: Count entries with underscore after modification
    after_with_underscore = df['Genus_CRBC_Updated'].astype(str).str.contains('_').sum()
    print(f"Entries with underscore after processing: {after_with_underscore}")
    
    # Step 5: Save output
    print(f"Saving results to: {args.output}")
    df.to_csv(args.output, index=False)
    
    print("\n=== Summary ===")
    print(f"Total rows processed: {len(df)}")
    print(f"Suffixes removed: {original_with_underscore - after_with_underscore}")
    print("Done!")


if __name__ == '__main__':
    main()

