#!/usr/bin/env python3
"""
Extract rows from bac120_metadata_r207.tsv where genus name matches those in Genus_Name.csv.
Filter for gtdb_representative='t' and keep only the first occurrence of each genus.

Usage:
    python 207_extract_matching_taxonomy.py -t METADATA_TSV -g GENUS_CSV -o OUTPUT

Example:
    python 207_extract_matching_taxonomy.py -t Result/NCBI_4395_Batch/05_Tree/Genus_Level/bac120_metadata_r207.tsv -g Result/NCBI_4395_Batch/05_Tree/Genus_Level/Genus_Name.csv -o Result/NCBI_4395_Batch/05_Tree/Genus_Level/filtered_metadata.tsv
"""

import argparse
import re
from tqdm import tqdm


def extract_genus_from_taxonomy(taxonomy_string):
    """
    Extract genus name from taxonomy string.
    
    Args:
        taxonomy_string: String like "d__Bacteria;p__...;g__Escherichia;s__..."
    
    Returns:
        Genus name (without g__ prefix) or None if not found
    """
    # Match g__ followed by any characters until ; or end of string
    match = re.search(r'g__([^;]+)', taxonomy_string)
    if match:
        return match.group(1)
    return None


def main():
    # Step 1: Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Extract rows from metadata TSV where gtdb_representative='t', genus matches target list, and deduplicate by genus."
    )
    parser.add_argument(
        "-t", "--taxonomy",
        required=True,
        help="Input metadata TSV file (bac120_metadata_r207.tsv)"
    )
    parser.add_argument(
        "-g", "--genus",
        required=True,
        help="Input genus CSV file (Genus_Name.csv)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file path"
    )
    
    args = parser.parse_args()
    
    # Step 2: Read target genus names from CSV
    print(f"Reading target genus names from: {args.genus}")
    target_genera = set()
    with open(args.genus, 'r') as f:
        for line in f:
            genus = line.strip()
            if genus:  # Skip empty lines
                target_genera.add(genus)
    print(f"Loaded {len(target_genera)} target genus names")
    
    # Step 3: Count total lines for progress bar
    print(f"Counting lines in metadata file...")
    with open(args.taxonomy, 'r') as f:
        total_lines = sum(1 for _ in f)
    print(f"Total lines: {total_lines}")
    
    # Step 4: Process metadata file and extract matching rows (deduplicated)
    print(f"Processing metadata file: {args.taxonomy}")
    matched_count = 0
    seen_genera = set()  # Track genera we've already written
    header_written = False
    
    with open(args.taxonomy, 'r') as f_in, open(args.output, 'w') as f_out:
        header = f_in.readline()
        f_out.write(header)  # Write header
        header_written = True
        
        # Find column indices
        header_parts = header.strip().split('\t')
        try:
            rep_idx = header_parts.index('gtdb_representative')
            tax_idx = header_parts.index('gtdb_taxonomy')
        except ValueError as e:
            print(f"Error: Required column not found - {e}")
            return
        
        for line in tqdm(f_in, total=total_lines-1, desc="Processing"):
            # Split line by tab to get columns
            parts = line.strip().split('\t')
            if len(parts) > max(rep_idx, tax_idx):
                gtdb_representative = parts[rep_idx]
                taxonomy = parts[tax_idx]
                
                # Check if gtdb_representative is 't'
                if gtdb_representative == 't':
                    # Extract genus name from taxonomy string
                    genus = extract_genus_from_taxonomy(taxonomy)
                    
                    # Check if genus matches target genera AND not yet written
                    if genus and genus in target_genera and genus not in seen_genera:
                        f_out.write(line)
                        seen_genera.add(genus)  # Mark this genus as written
                        matched_count += 1
    
    # Step 5: Print summary
    print(f"\nResults saved to: {args.output}")
    print(f"Total unique genera matched: {matched_count}")
    print(f"Target genera count: {len(target_genera)}")
    if matched_count < len(target_genera):
        missing_count = len(target_genera) - matched_count
        print(f"Warning: {missing_count} genera from target list not found in metadata file")


if __name__ == "__main__":
    main()

