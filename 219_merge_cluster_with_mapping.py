#!/usr/bin/env python3
"""
Merge cluster rename map with contig sample mapping data.
Extract old contig ID from original column and match with Contig_ID in mapping file.
"""

import argparse
import pandas as pd
from tqdm import tqdm


def count_values(value):
    """Count number of values separated by commas. Empty value counts as 0."""
    if pd.isna(value) or value == '':
        return 0
    return value.count(',') + 1


def main():
    parser = argparse.ArgumentParser(
        description='Merge cluster rename map with contig sample mapping data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 219_merge_cluster_with_mapping.py -i Result/NCBI_4395_Batch/07_Network/coocc_network/mmseq_overall_contigs_cluster.rename_map.tsv -m Result/NCBI_4395_Batch/07_Network/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv -o Result/NCBI_4395_Batch/07_Network/coocc_network/mmseq_overall_contigs_cluster.rename_map.merged.tsv
        """
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input TSV file (mmseq_overall_contigs_cluster.rename_map.tsv)')
    parser.add_argument('-m', '--mapping', required=True,
                        help='Mapping CSV file (Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output TSV file')
    
    args = parser.parse_args()
    
    # Step 1: Read cluster rename map
    print("Reading cluster rename map...")
    df_cluster = pd.read_csv(args.input, sep='\t')
    
    # Step 2: Extract old contig ID from original column
    # Remove content before first | (including the |)
    print("Extracting old contig IDs...")
    df_cluster['old'] = df_cluster['original'].str.split('|', n=1).str[1]
    
    # Step 3: Read mapping file
    print("Reading mapping file...")
    df_mapping = pd.read_csv(args.mapping, low_memory=False)
    
    # Step 4: Select required columns from mapping file
    required_cols = [
        'Contig_ID', 'Sample_ID', 'Contig_Type', 'Defense_Subtype',
        'AntiDS_Type', 'AMR_Type', 'Host', 'Kingdom_CRBC', 'Phylum_CRBC',
        'Class_CRBC', 'Order_CRBC', 'Family_CRBC', 'Genus_CRBC',
        'Species_CRBC', 'Genus_CRBC_Updated', 'predicted_mobility'
    ]
    
    # Check if all required columns exist
    missing_cols = [col for col in required_cols if col not in df_mapping.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in mapping file: {missing_cols}")
    
    df_mapping_subset = df_mapping[required_cols].copy()
    
    # Step 5: Merge based on old contig ID and Contig_ID
    print("Merging data...")
    df_merged = df_cluster.merge(
        df_mapping_subset,
        left_on='old',
        right_on='Contig_ID',
        how='left'
    )
    
    # Step 6: Calculate count columns
    print("Calculating count columns...")
    df_merged['Defense_Num'] = df_merged['Defense_Subtype'].apply(count_values)
    df_merged['AntiDS_Num'] = df_merged['AntiDS_Type'].apply(count_values)
    df_merged['AMR_Num'] = df_merged['AMR_Type'].apply(count_values)
    
    # Step 7: Write output
    print(f"Writing output to {args.output}...")
    df_merged.to_csv(args.output, sep='\t', index=False)
    
    print(f"Done! Merged {len(df_merged)} rows.")
    print(f"Matched rows: {df_merged['Sample_ID'].notna().sum()}")
    print(f"Unmatched rows: {df_merged['Sample_ID'].isna().sum()}")


if __name__ == '__main__':
    main()

