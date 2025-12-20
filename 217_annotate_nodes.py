#!/usr/bin/env python3
"""
Annotate nodes files with information from merged cluster rename map.
Match new column from merged file with cluster/contig column from nodes files.
"""

import argparse
import pandas as pd
import os
from pathlib import Path
from tqdm import tqdm


# Color mapping based on Class_CRBC
CLASS_COLOR_MAP = {
    'Actinomycetia': '#98df8a',
    'Alphaproteobacteria': '#aec7e8',
    'Bacilli': '#ff7f0e',
    'Thermoleophilia': '#ff9896',
    'Bacteroidia': '#d62728',
    'Gammaproteobacteria': '#ffbb78',
    'Deinococci': '#1f77b4',
    'Acidimicrobiia': '#2ca02c',
    'Campylobacteria': '#9467bd'
}


def main():
    parser = argparse.ArgumentParser(
        description='Annotate nodes files with cluster mapping data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # For coocc_network (using cluster column)
  python 217_annotate_nodes.py -m Result/NCBI_4395_Batch/07_Network/coocc_network/mmseq_overall_contigs_cluster.rename_map.merged.tsv -d Result/NCBI_4395_Batch/07_Network/coocc_network -o Result/NCBI_4395_Batch/07_Network/coocc_network/Annotation -k cluster
  
  # For transfer_network (using contig column)
  python 217_annotate_nodes.py -m Result/NCBI_4395_Batch/07_Network/transfer_network/mmseq_overall_contigs_cluster.rename_map.merged.tsv -d Result/NCBI_4395_Batch/07_Network/transfer_network -o Result/NCBI_4395_Batch/07_Network/transfer_network/Annotation -k contig
        """
    )
    parser.add_argument('-m', '--merged', required=True,
                        help='Merged cluster rename map TSV file')
    parser.add_argument('-d', '--directory', required=True,
                        help='Directory containing nodes TSV files')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Output directory (Annotation folder)')
    parser.add_argument('-k', '--key-column', default='auto',
                        choices=['auto', 'cluster', 'contig'],
                        help='Column name to use for matching (auto: detect from nodes file, cluster: use cluster column, contig: use contig column)')
    
    args = parser.parse_args()
    
    # Step 1: Read merged file
    print(f"Reading merged file: {args.merged}...")
    df_merged = pd.read_csv(args.merged, sep='\t', low_memory=False)
    
    # Step 2: Select columns to add
    columns_to_add = [
        'new', 'Sample_ID', 'Contig_Type', 'Defense_Subtype',
        'AntiDS_Type', 'AMR_Type', 'Host', 'Kingdom_CRBC', 'Phylum_CRBC',
        'Class_CRBC', 'Order_CRBC', 'Family_CRBC', 'Genus_CRBC',
        'Species_CRBC', 'Genus_CRBC_Updated', 'predicted_mobility',
        'Defense_Num', 'AntiDS_Num', 'AMR_Num'
    ]
    
    # Check if all columns exist
    missing_cols = [col for col in columns_to_add if col not in df_merged.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in merged file: {missing_cols}")
    
    df_merged_subset = df_merged[columns_to_add].copy()
    
    # Step 3: Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 4: Find all nodes TSV files
    nodes_files = list(Path(args.directory).glob('*_nodes.tsv'))
    
    if not nodes_files:
        print(f"No nodes TSV files found in {args.directory}")
        return
    
    print(f"Found {len(nodes_files)} nodes file(s)")
    
    # Step 5: Process each nodes file
    for nodes_file in tqdm(nodes_files, desc="Processing nodes files"):
        print(f"\nProcessing: {nodes_file.name}")
        
        # Read nodes file
        df_nodes = pd.read_csv(nodes_file, sep='\t', low_memory=False)
        
        # Determine key column for matching
        key_col = args.key_column
        if key_col == 'auto':
            # Auto-detect: prefer cluster, then contig
            if 'cluster' in df_nodes.columns:
                key_col = 'cluster'
            elif 'contig' in df_nodes.columns:
                key_col = 'contig'
            else:
                print(f"  Warning: Neither 'cluster' nor 'contig' column found in {nodes_file.name}, skipping...")
                continue
        
        # Check if key column exists
        if key_col not in df_nodes.columns:
            print(f"  Warning: '{key_col}' column not found in {nodes_file.name}, skipping...")
            continue
        
        # Merge based on new (from merged) and key_col (from nodes)
        df_annotated = df_nodes.merge(
            df_merged_subset,
            left_on=key_col,
            right_on='new',
            how='left'
        )
        
        # Drop the 'new' column as it's redundant with key_col
        if 'new' in df_annotated.columns:
            df_annotated = df_annotated.drop(columns=['new'])
        
        # Rename key_col to id for output
        if key_col in df_annotated.columns:
            df_annotated = df_annotated.rename(columns={key_col: 'id'})
        
        # Add color column based on Class_CRBC
        if 'Class_CRBC' in df_annotated.columns:
            df_annotated['color'] = df_annotated['Class_CRBC'].map(CLASS_COLOR_MAP)
        else:
            df_annotated['color'] = None
        
        # Rename Defense_Subtype to label
        if 'Defense_Subtype' in df_annotated.columns:
            df_annotated = df_annotated.rename(columns={'Defense_Subtype': 'label'})
        
        # Save to output directory
        output_file = output_dir / nodes_file.name
        df_annotated.to_csv(output_file, sep='\t', index=False)
        
        matched_count = df_annotated['Sample_ID'].notna().sum()
        total_count = len(df_annotated)
        print(f"  Saved to: {output_file}")
        print(f"  Matched: {matched_count}/{total_count} rows")
    
    print("\nDone!")


if __name__ == '__main__':
    main()

