#!/usr/bin/env python3
"""
Filter Sample_Genus_map.csv based on Genus_CRBC_Updated values in intersection file.
Extract rows from Sample_Genus_map.csv where Genus_CRBC_Updated matches values in intersection file.

Usage:
    python 209_filter_sample_genus_map.py -i INTERSECTION -s SAMPLE_MAP -o OUTPUT

Example:
    python 209_filter_sample_genus_map.py -i Result/NCBI_4395_Batch/05_Tree/Genus_Level_Intersection/Genus_CRBC_Updated_Intersection.csv -s Result/NCBI_4395_Batch/05_Tree/Genus_Level_Intersection/Sample_Genus_map.csv -o Result/NCBI_4395_Batch/05_Tree/Genus_Level_Intersection/Sample_Genus_map_Intersection.csv
"""

import argparse
import pandas as pd
from tqdm import tqdm

def filter_sample_genus_map(intersection_file, sample_genus_map_file, output_file):
    """
    Filter Sample_Genus_map.csv based on Genus_CRBC_Updated values in intersection file.
    
    Args:
        intersection_file: Input CSV file with Genus_CRBC_Updated values to match
        sample_genus_map_file: Input CSV file (Sample_Genus_map.csv) to filter
        output_file: Output CSV file path
    """
    # Step 1: Read intersection file to get Genus_CRBC_Updated values
    print(f"Reading intersection file: {intersection_file}")
    intersection_df = pd.read_csv(intersection_file)
    intersection_genus = set(intersection_df['Genus_CRBC_Updated'].dropna().unique())
    print(f"Found {len(intersection_genus)} unique Genus_CRBC_Updated values in intersection")
    
    # Step 2: Read Sample_Genus_map file
    print(f"\nReading sample-genus map file: {sample_genus_map_file}")
    sample_genus_df = pd.read_csv(sample_genus_map_file)
    print(f"Total rows in sample-genus map: {len(sample_genus_df)}")
    
    # Step 3: Filter rows where Genus_CRBC_Updated is in intersection set
    print("\nFiltering rows...")
    filtered_df = sample_genus_df[sample_genus_df['Genus_CRBC_Updated'].isin(intersection_genus)].copy()
    print(f"Found {len(filtered_df)} matching rows")
    
    # Step 4: Sort by Genus_CRBC_Updated, then by Sample_ID
    filtered_df = filtered_df.sort_values(['Genus_CRBC_Updated', 'Sample_ID']).reset_index(drop=True)
    
    # Step 5: Save to output file
    print(f"\nSaving filtered results to: {output_file}")
    filtered_df.to_csv(output_file, index=False)
    
    print(f"Done! Saved {len(filtered_df)} entries to {output_file}")
    
    # Step 6: Print summary statistics
    print("\n=== Summary ===")
    print(f"Unique Genus_CRBC_Updated in intersection: {len(intersection_genus)}")
    print(f"Total rows in original sample-genus map: {len(sample_genus_df)}")
    print(f"Filtered rows: {len(filtered_df)}")
    
    # Count samples per genus
    if len(filtered_df) > 0:
        genus_counts = filtered_df['Genus_CRBC_Updated'].value_counts().sort_index()
        print(f"\nSamples per Genus_CRBC_Updated:")
        for genus, count in genus_counts.items():
            print(f"  {genus}: {count}")

def main():
    parser = argparse.ArgumentParser(
        description='Filter Sample_Genus_map.csv based on Genus_CRBC_Updated values in intersection file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 209_filter_sample_genus_map.py -i intersection.csv -s sample_genus_map.csv -o filtered_output.csv
  python 209_filter_sample_genus_map.py --intersection intersection.csv --sample-map sample_genus_map.csv --output filtered_output.csv
        """
    )
    parser.add_argument('-i', '--intersection', 
                        required=True,
                        help='Input CSV file with Genus_CRBC_Updated intersection values')
    parser.add_argument('-s', '--sample-map',
                        required=True,
                        help='Input CSV file (Sample_Genus_map.csv) to filter')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output CSV file path')
    
    args = parser.parse_args()
    filter_sample_genus_map(args.intersection, args.sample_map, args.output)

if __name__ == '__main__':
    main()

