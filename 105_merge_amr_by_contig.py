#!/usr/bin/env python3
"""
Merge AMR results by Contig_ID.
For each Contig_ID, merge all entries and concatenate values with commas (including duplicates).
"""

import argparse
import pandas as pd
from tqdm import tqdm

def merge_amr_by_contig(input_file, output_file):
    """
    Merge AMR results by Contig_ID.
    
    Args:
        input_file: Input CSV file path
        output_file: Output CSV file path
    """
    # Step 1: Read input CSV file
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file)
    
    print(f"Total rows: {len(df)}")
    print(f"Unique Contig_IDs: {df['Contig_ID'].nunique()}")
    
    # Step 2: Group by Contig_ID and merge
    print("Merging by Contig_ID...")
    
    # Get all column names except Contig_ID
    merge_columns = [col for col in df.columns if col != 'Contig_ID']
    
    # Define aggregation functions for each column
    def make_agg_func(col_name):
        if col_name in ['AMR_Start', 'AMR_End']:
            # Keep as float for numeric columns
            return lambda x: ','.join(x.astype(float).astype(str))
        else:
            # Join with comma for other columns
            return lambda x: ','.join(map(str, x))
    
    agg_dict = {col: make_agg_func(col) for col in merge_columns}
    
    # Group by Contig_ID and aggregate
    merged = df.groupby('Contig_ID', as_index=False).agg(agg_dict)
    
    print(f"Merged rows: {len(merged)}")
    
    # Step 3: Save to output file
    print(f"Saving to: {output_file}")
    merged.to_csv(output_file, index=False)
    
    print("Done!")

def main():
    parser = argparse.ArgumentParser(
        description='Merge AMR results by Contig_ID, concatenating values with commas (including duplicates).'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CSV file with AMR results'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output CSV file with merged results'
    )
    
    args = parser.parse_args()
    
    merge_amr_by_contig(args.input, args.output)

if __name__ == '__main__':
    main()

