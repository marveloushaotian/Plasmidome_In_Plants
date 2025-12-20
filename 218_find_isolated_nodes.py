#!/usr/bin/env python3
"""
Find isolated nodes (nodes that appear in nodes files but not in edges files).
These are black nodes with no connections.
"""

import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm


def find_isolated_nodes(nodes_file, edges_file, node_id_col=None, output_file=None):
    """
    Find nodes that appear in nodes file but not in edges file.
    
    Args:
        nodes_file: Path to nodes TSV file
        edges_file: Path to edges TSV file
        node_id_col: Column name for node ID in nodes file (auto-detect if None)
        output_file: Optional path to save isolated nodes
    
    Returns:
        DataFrame with isolated nodes
    """
    # Read nodes file
    df_nodes = pd.read_csv(nodes_file, sep='\t', low_memory=False)
    
    # Auto-detect node ID column if not specified
    if node_id_col is None:
        possible_cols = ['id', 'cluster', 'contig', 'node']
        for col in possible_cols:
            if col in df_nodes.columns:
                node_id_col = col
                break
        if node_id_col is None:
            raise ValueError(f"Could not auto-detect node ID column. Available columns: {list(df_nodes.columns)}")
    
    if node_id_col not in df_nodes.columns:
        raise ValueError(f"Column '{node_id_col}' not found in {nodes_file}. Available columns: {list(df_nodes.columns)}")
    
    # Get all node IDs from nodes file
    all_node_ids = set(df_nodes[node_id_col].unique())
    
    # Read edges file
    df_edges = pd.read_csv(edges_file, sep='\t', low_memory=False)
    
    # Get all node IDs that appear in edges (from source and target columns)
    connected_node_ids = set()
    if 'source' in df_edges.columns:
        connected_node_ids.update(df_edges['source'].unique())
    if 'target' in df_edges.columns:
        connected_node_ids.update(df_edges['target'].unique())
    
    # Find isolated nodes (in nodes but not in edges)
    isolated_node_ids = all_node_ids - connected_node_ids
    
    # Get full information for isolated nodes
    df_isolated = df_nodes[df_nodes[node_id_col].isin(isolated_node_ids)].copy()
    
    # Save if output file is specified
    if output_file:
        df_isolated.to_csv(output_file, sep='\t', index=False)
    
    return df_isolated, len(isolated_node_ids), len(all_node_ids)


def main():
    parser = argparse.ArgumentParser(
        description='Find isolated nodes in network files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # For coocc_network
  python 218_find_isolated_nodes.py -d Result/NCBI_4395_Batch/07_Network/coocc_network/Annotation -n cluster -o Result/NCBI_4395_Batch/07_Network/coocc_network/isolated_nodes
  
  # For transfer_network
  python 218_find_isolated_nodes.py -d Result/NCBI_4395_Batch/07_Network/transfer_network/Annotation -n contig -o Result/NCBI_4395_Batch/07_Network/transfer_network/isolated_nodes
        """
    )
    parser.add_argument('-d', '--directory', required=True,
                        help='Directory containing nodes and edges TSV files')
    parser.add_argument('-n', '--node-id-col', default=None,
                        help='Column name for node ID in nodes files (auto-detect if not specified: id, cluster, contig, node)')
    parser.add_argument('-o', '--output-dir', default='isolated_nodes',
                        help='Output directory for isolated nodes files (default: isolated_nodes)')
    
    args = parser.parse_args()
    
    # Find all nodes files
    nodes_files = sorted(Path(args.directory).glob('*_nodes.tsv'))
    
    if not nodes_files:
        print(f"No nodes TSV files found in {args.directory}")
        return
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Found {len(nodes_files)} nodes file(s)\n")
    
    # Process each nodes file
    results = []
    for nodes_file in tqdm(nodes_files, desc="Processing files"):
        # Find corresponding edges file
        edges_file = nodes_file.parent / nodes_file.name.replace('_nodes.tsv', '_edges.tsv')
        
        if not edges_file.exists():
            print(f"  Warning: Edges file not found for {nodes_file.name}, skipping...")
            continue
        
        print(f"\nProcessing: {nodes_file.name}")
        
        # Find isolated nodes
        output_file = output_dir / nodes_file.name.replace('_nodes.tsv', '_isolated_nodes.tsv')
        try:
            # Try to auto-detect node ID column for this file
            df_test = pd.read_csv(nodes_file, sep='\t', nrows=1)
            possible_cols = ['id', 'cluster', 'contig', 'node']
            node_id_col = args.node_id_col
            if node_id_col is None:
                for col in possible_cols:
                    if col in df_test.columns:
                        node_id_col = col
                        break
            
            if node_id_col is None:
                print(f"  Error: Could not auto-detect node ID column in {nodes_file.name}")
                continue
            
            df_isolated, isolated_count, total_count = find_isolated_nodes(
                nodes_file, edges_file, node_id_col, output_file
            )
            
            percentage = (isolated_count / total_count * 100) if total_count > 0 else 0
            results.append({
                'File': nodes_file.name,
                'Total_Nodes': total_count,
                'Isolated_Nodes': isolated_count,
                'Connected_Nodes': total_count - isolated_count,
                'Percentage_Isolated': f"{percentage:.2f}%"
            })
            
            print(f"  Total nodes: {total_count}")
            print(f"  Isolated nodes: {isolated_count} ({percentage:.2f}%)")
            print(f"  Connected nodes: {total_count - isolated_count}")
            print(f"  Saved to: {output_file}")
            
        except Exception as e:
            print(f"  Error processing {nodes_file.name}: {e}")
    
    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    df_summary = pd.DataFrame(results)
    if not df_summary.empty:
        print(df_summary.to_string(index=False))
        
        # Save summary
        summary_file = output_dir / 'isolated_nodes_summary.tsv'
        df_summary.to_csv(summary_file, sep='\t', index=False)
        print(f"\nSummary saved to: {summary_file}")
    
    print("\nDone!")


if __name__ == '__main__':
    main()

