#!/usr/bin/env python3
"""
Extract intersection of Genus_CRBC_Updated values across different Host groups.
For each unique Host, get all unique Genus_CRBC_Updated values,
then calculate the intersection (values that appear in all Host groups).

Usage:
    python 210_get_genus_intersection.py -i INPUT -o OUTPUT

Example:
    python 210_get_genus_intersection.py -i Result/NCBI_4395_Batch/05_Tree/Genus_Level_Intersection/Sample_Genus_map.csv -o Result/NCBI_4395_Batch/05_Tree/Genus_Level_Intersection/Genus_CRBC_Updated_Intersection.csv
"""

import argparse
import pandas as pd
from tqdm import tqdm

def get_genus_intersection(input_file, output_file):
    """
    Get intersection of Genus_CRBC_Updated values across all Host groups.
    
    Args:
        input_file: Input CSV file path
        output_file: Output CSV file path
    """
    # Step 1: Read input CSV file
    print(f"Reading input file: {input_file}")
    df = pd.read_csv(input_file)
    
    print(f"Total rows: {len(df)}")
    
    # Step 2: Get unique Host values
    unique_hosts = df['Host'].dropna().unique()
    print(f"Found {len(unique_hosts)} unique Host values: {sorted(unique_hosts)}")
    
    # Step 3: Get unique Genus_CRBC_Updated for each Host
    host_genus_dict = {}
    for host in tqdm(unique_hosts, desc="Processing Host groups"):
        host_df = df[df['Host'] == host]
        # Get unique Genus_CRBC_Updated values (drop NaN and empty strings)
        genus_values = host_df['Genus_CRBC_Updated'].dropna()
        genus_values = genus_values[genus_values != ''].unique()
        host_genus_dict[host] = set(genus_values)
        print(f"  {host}: {len(genus_values)} unique Genus_CRBC_Updated values")
    
    # Step 4: Calculate intersection across all Host groups
    print("\nCalculating intersection...")
    if not host_genus_dict:
        print("Warning: No Host groups found!")
        intersection = set()
    else:
        # Start with first host's set, then intersect with others
        intersection = host_genus_dict[list(host_genus_dict.keys())[0]]
        for host in list(host_genus_dict.keys())[1:]:
            intersection = intersection & host_genus_dict[host]
    
    print(f"Intersection contains {len(intersection)} unique Genus_CRBC_Updated values")
    
    # Step 5: Get taxonomy information for each Genus_CRBC_Updated in intersection
    print("\nExtracting taxonomy information...")
    taxonomy_cols = ['Kingdom_CRBC', 'Phylum_CRBC', 'Class_CRBC', 'Order_CRBC', 'Family_CRBC', 'Genus_CRBC_Updated']
    
    if intersection:
        # Filter dataframe for Genus_CRBC_Updated values in intersection
        intersection_df = df[df['Genus_CRBC_Updated'].isin(intersection) & 
                            df['Genus_CRBC_Updated'].notna() & 
                            (df['Genus_CRBC_Updated'] != '')]
        
        # Get unique combinations (for each Genus_CRBC_Updated, get its taxonomy info)
        # Group by Genus_CRBC_Updated and take first occurrence of each
        result_df = intersection_df[taxonomy_cols].drop_duplicates(subset=['Genus_CRBC_Updated'])
        
        # Sort by Genus_CRBC_Updated
        result_df = result_df.sort_values('Genus_CRBC_Updated').reset_index(drop=True)
        
        # Reorder columns: Genus_CRBC_Updated first, then other taxonomy columns
        result_df = result_df[['Genus_CRBC_Updated', 'Kingdom_CRBC', 'Phylum_CRBC', 'Class_CRBC', 'Order_CRBC', 'Family_CRBC']]
    else:
        result_df = pd.DataFrame(columns=['Genus_CRBC_Updated', 'Kingdom_CRBC', 'Phylum_CRBC', 'Class_CRBC', 'Order_CRBC', 'Family_CRBC'])
    
    # Step 6: Save to output file
    print(f"\nSaving results to: {output_file}")
    result_df.to_csv(output_file, index=False)
    
    print(f"Done! Saved {len(result_df)} entries to {output_file}")
    
    # Step 7: Print summary statistics
    print("\n=== Summary ===")
    print(f"Total Host groups: {len(unique_hosts)}")
    print(f"Genus_CRBC_Updated values in intersection: {len(intersection)}")
    if intersection:
        print(f"\nIntersection values:")
        for genus in sorted(intersection):
            print(f"  - {genus}")

def main():
    parser = argparse.ArgumentParser(
        description='Extract intersection of Genus_CRBC_Updated values across all Host groups',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 210_get_genus_intersection.py -i input.csv -o output.csv
  python 210_get_genus_intersection.py --input data.csv --output intersection.csv
        """
    )
    parser.add_argument('-i', '--input', 
                        required=True,
                        help='Input CSV file path')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='Output CSV file path')
    
    args = parser.parse_args()
    get_genus_intersection(args.input, args.output)

if __name__ == '__main__':
    main()

