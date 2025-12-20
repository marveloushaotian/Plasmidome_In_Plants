#!/usr/bin/env python3
"""
Prune phylogenetic tree and restore internal node labels.

This script:
1. Prunes the GTDB tree to keep only target genera
2. Renames leaf nodes to genus names only
3. Restores internal node labels (phylum, class, order, family) based on taxonomy

This ensures that:
- Leaf nodes show only genus names (e.g., "Escherichia")
- Internal nodes show taxonomic labels (e.g., "p__Proteobacteria", "f__Enterobacteriaceae")
- Branch lengths (evolutionary distances) are preserved

Usage:
    python 211_prune_and_restore_labels.py -t TREE_FILE -i ID_FILE -x TAXONOMY_FILE -o OUTPUT

Example:
    python 211_prune_and_restore_labels.py -t Result/NCBI_4395_Batch/05_Tree/Genus_Level/bac120_r207.tree -i Result/NCBI_4395_Batch/05_Tree/Genus_Level/Archive/filtered_bac120_metadata_r207_unique_id.tsv -x Result/NCBI_4395_Batch/05_Tree/Genus_Level/Archive/filtered_bac120_metadata_r207_unique_id_alllevel.csv -o Result/NCBI_4395_Batch/05_Tree/Genus_Level/Archive/pruned_tree_with_internal_labels.tree
"""

import csv
import argparse
from collections import defaultdict
from ete3 import Tree


def read_id_list(id_file):
    """Read genome IDs from TSV file."""
    ids = []
    with open(id_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                ids.append(line)
    return ids


def read_full_taxonomy(taxonomy_file):
    """
    Read full taxonomy for each genome ID.
    
    Returns:
        dict: {genome_id: {
            'genus': genus_name,
            'domain': d__...,
            'phylum': p__...,
            'class': c__...,
            'order': o__...,
            'family': f__...
        }}
    """
    taxonomy = {}
    with open(taxonomy_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row['Acc'].strip()
            taxonomy[acc] = {
                'genus': row['GTDB_Genus'].strip().replace('g__', ''),  # Remove g__ prefix
                'domain': row['GTDB_Domain'].strip(),
                'phylum': row['GTDB_Phylum'].strip(),
                'class': row['GTDB_Class'].strip(),
                'order': row['GTDB_Order'].strip(),
                'family': row['GTDB_Family'].strip()
            }
    return taxonomy


def find_monophyletic_groups(tree, genome_taxonomy):
    """
    Find monophyletic groups in the tree and determine which internal nodes
    should have which taxonomic labels.
    
    Returns:
        dict: {node: label} mapping for internal nodes
    """
    node_labels = {}
    
    # For each internal node, check if all its descendants belong to the same taxon
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            continue
        
        # Get all leaf descendants
        leaves = node.get_leaves()
        leaf_names = [leaf.name for leaf in leaves]
        
        # Check taxonomy for each taxonomic level
        for rank in ['phylum', 'class', 'order', 'family']:
            taxa_at_rank = set()
            for leaf_name in leaf_names:
                if leaf_name in genome_taxonomy:
                    taxa_at_rank.add(genome_taxonomy[leaf_name][rank])
            
            # If all leaves belong to the same taxon at this rank
            if len(taxa_at_rank) == 1:
                taxon = taxa_at_rank.pop()
                if taxon:  # Not empty
                    # This node represents this taxon
                    # Store with the highest resolution label
                    if rank == 'family' or node not in node_labels:
                        node_labels[node] = taxon
    
    return node_labels


def prune_and_restore(tree_file, id_file, taxonomy_file, output_file):
    """
    Main function to prune tree and restore internal labels.
    """
    print("=" * 80)
    print("   Pruning Tree and Restoring Internal Node Labels")
    print("=" * 80)
    print()
    
    # Step 1: Read genome IDs
    print("Step 1: Reading genome IDs...")
    target_ids = read_id_list(id_file)
    print(f"  Loaded {len(target_ids)} genome IDs")
    print()
    
    # Step 2: Read taxonomy
    print("Step 2: Reading taxonomy information...")
    genome_taxonomy = read_full_taxonomy(taxonomy_file)
    print(f"  Loaded taxonomy for {len(genome_taxonomy)} genomes")
    print()
    
    # Step 3: Load and prune tree
    print("Step 3: Loading original tree...")
    tree = Tree(tree_file, format=1, quoted_node_names=True)
    print(f"  Original tree: {len(tree.get_leaves())} leaves")
    
    # Find which target IDs are in the tree
    tree_leaves = set([leaf.name for leaf in tree.get_leaves()])
    present_ids = [gid for gid in target_ids if gid in tree_leaves]
    
    print(f"  Target genomes in tree: {len(present_ids)}")
    print()
    
    print("Step 4: Pruning tree...")
    tree.prune(present_ids, preserve_branch_length=True)
    print(f"  Pruned tree: {len(tree.get_leaves())} leaves")
    print()
    
    # Step 5: Rename leaves to genus names only
    print("Step 5: Renaming leaves to genus names...")
    renamed_count = 0
    # Create a mapping from genome ID to taxonomy for remaining leaves
    leaf_taxonomy = {}
    
    for leaf in tree.get_leaves():
        genome_id = leaf.name
        if genome_id in genome_taxonomy:
            genus_name = genome_taxonomy[genome_id]['genus']
            leaf.name = genus_name
            leaf_taxonomy[genus_name] = genome_taxonomy[genome_id]
            renamed_count += 1
    
    print(f"  Renamed {renamed_count} leaves to genus names")
    print()
    
    # Step 6: Find and restore internal node labels
    print("Step 6: Restoring internal node labels...")
    
    # Re-create taxonomy mapping with new leaf names
    leaf_name_taxonomy = {}
    for leaf in tree.get_leaves():
        genus_name = leaf.name
        if genus_name in leaf_taxonomy:
            leaf_name_taxonomy[genus_name] = leaf_taxonomy[genus_name]
    
    # Find monophyletic groups
    node_labels = find_monophyletic_groups(tree, leaf_name_taxonomy)
    
    # Apply labels to nodes
    for node, label in node_labels.items():
        node.name = label
    
    print(f"  Restored {len(node_labels)} internal node labels")
    print()
    
    # Count labels by rank
    rank_counts = defaultdict(int)
    for label in node_labels.values():
        if label.startswith('p__'):
            rank_counts['phylum'] += 1
        elif label.startswith('c__'):
            rank_counts['class'] += 1
        elif label.startswith('o__'):
            rank_counts['order'] += 1
        elif label.startswith('f__'):
            rank_counts['family'] += 1
    
    print("  Labels by rank:")
    print(f"    Phylum (p__):  {rank_counts['phylum']}")
    print(f"    Class (c__):   {rank_counts['class']}")
    print(f"    Order (o__):   {rank_counts['order']}")
    print(f"    Family (f__):  {rank_counts['family']}")
    print()
    
    # Step 7: Write tree
    print("Step 7: Writing tree...")
    tree.write(outfile=output_file, format=1)
    print(f"  Wrote: {output_file}")
    print()
    
    print("=" * 80)
    print("Success!")
    print()
    print("The tree now contains:")
    print(f"  • {len(tree.get_leaves())} genera as leaf nodes (genus names only)")
    print(f"  • {len(node_labels)} internal nodes with taxonomic labels")
    print("  • Preserved branch lengths (evolutionary distances)")
    print()
    print("Leaf node format: Escherichia, Bacillus, Pseudomonas, ...")
    print("Internal node format: p__Proteobacteria, c__Gammaproteobacteria, ...")
    print("=" * 80)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Prune phylogenetic tree and restore internal node labels",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '-t', '--tree',
        required=True,
        help='Input tree file (bac120_r207.tree)'
    )
    parser.add_argument(
        '-i', '--id-file',
        required=True,
        help='Input ID file (filtered_bac120_metadata_r207_unique_id.tsv)'
    )
    parser.add_argument(
        '-x', '--taxonomy-file',
        required=True,
        help='Input taxonomy file (filtered_bac120_metadata_r207_unique_id_alllevel.csv)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output tree file'
    )
    
    args = parser.parse_args()
    
    prune_and_restore(args.tree, args.id_file, args.taxonomy_file, args.output)


if __name__ == "__main__":
    main()

