#!/usr/bin/env python3
"""
Filter rows from Contig_Sample_Mapping file based on Genus names.

This script extracts rows where Genus_CRBC_Updated matches genera in Genus_Name.csv,
and outputs only Sample_ID and Genus_CRBC_Updated columns.

Usage:
    python 208_filter_genus_samples.py -i INPUT -g GENUS -o OUTPUT

Example:
    python 208_filter_genus_samples.py -i Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv -g Result/NCBI_4395_Batch/05_Tree/Genus_Level/Genus_Name.csv -o Result/NCBI_4395_Batch/05_Tree/Genus_Level/Filtered_Sample_Genus.csv
"""

import pandas as pd
import argparse
from tqdm import tqdm

def main():
    # Step 1: Parse arguments
    parser = argparse.ArgumentParser(
        description='Filter Contig_Sample_Mapping by Genus names and keep only Sample_ID and Genus_CRBC_Updated columns.'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CSV file with contig sample mapping'
    )
    parser.add_argument(
        '-g', '--genus',
        required=True,
        help='Input CSV file with genus names'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output CSV file'
    )
    
    args = parser.parse_args()
    
    # Step 2: Read genus names file
    print(f"Reading genus names from {args.genus}...")
    genus_df = pd.read_csv(args.genus, header=None, names=['Genus'])
    genus_list = genus_df['Genus'].unique().tolist()
    print(f"Found {len(genus_list)} unique genera")
    
    # Step 3: Read main mapping file
    print(f"Reading mapping file from {args.input}...")
    mapping_df = pd.read_csv(args.input)
    print(f"Total rows in mapping file: {len(mapping_df)}")
    
    # Step 4: Filter rows where Genus_CRBC_Updated is in genus_list
    print("Filtering rows based on genus names...")
    filtered_df = mapping_df[mapping_df['Genus_CRBC_Updated'].isin(genus_list)]
    print(f"Rows after filtering: {len(filtered_df)}")
    
    # Step 5: Keep only Sample_ID and Genus_CRBC_Updated columns
    result_df = filtered_df[['Sample_ID', 'Genus_CRBC_Updated']].copy()
    
    # Step 6: Remove duplicates
    print("Removing duplicate rows...")
    print(f"Rows before deduplication: {len(result_df)}")
    result_df = result_df.drop_duplicates()
    print(f"Rows after deduplication: {len(result_df)}")
    
    # Step 7: Save output
    print(f"Writing output to {args.output}...")
    result_df.to_csv(args.output, index=False)
    print(f"Done! Output saved to {args.output}")
    print(f"Output contains {len(result_df)} rows and {len(result_df.columns)} columns")

if __name__ == '__main__':
    main()

