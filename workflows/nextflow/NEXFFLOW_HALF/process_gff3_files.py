#!/usr/bin/env python3
"""
Process GFF3 files to extract protein information.
Extracts sequence ID, protein ID, start and end positions.
"""

import os
import re
import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm

def extract_id_from_attributes(attributes):
    """
    Extract ID from the attributes column (column 9).
    
    Args:
        attributes: String like 'ID=1_1;partial=00;...'
    
    Returns:
        Last part of ID after underscore (e.g., '1' from '1_1')
    """
    # Step 1: Extract ID value using regex
    match = re.search(r'ID=([^;]+)', attributes)
    if match:
        id_value = match.group(1)
        # Step 2: Extract only the part after the last underscore
        parts = id_value.split('_')
        if len(parts) > 1:
            return parts[-1]  # Return last part after underscore
        return id_value
    return None

def process_gff3_file(gff3_file):
    """
    Process a single GFF3 file and extract relevant information.
    
    Args:
        gff3_file: Path to GFF3 file
    
    Returns:
        List of dictionaries with Prot_ID, Prot_Start, Prot_End
    """
    results = []
    
    with open(gff3_file, 'r') as f:
        for line in f:
            # Step 2: Skip comment lines (lines starting with #)
            if line.startswith('#'):
                continue
            
            # Step 3: Parse non-comment lines
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            # Step 4: Extract required columns
            seq_id = parts[0]          # Column 1: sequence ID
            start = parts[3]            # Column 4: start position
            end = parts[4]              # Column 5: end position
            attributes = parts[8]       # Column 9: attributes
            
            # Step 5: Extract ID from attributes
            prot_id_suffix = extract_id_from_attributes(attributes)
            if prot_id_suffix:
                # Step 6: Combine sequence ID with protein ID
                prot_id = f"{seq_id}_{prot_id_suffix}"
                
                results.append({
                    'Prot_ID': prot_id,
                    'Prot_Start': start,
                    'Prot_End': end
                })
    
    return results

def process_all_gff3_files(input_dir, output_file):
    """
    Process all GFF3 files in the directory.
    
    Args:
        input_dir: Directory containing GFF3 files
        output_file: Path to output CSV file
    """
    # Step 1: Find all GFF3 files
    print("Searching for .gff3 files...")
    gff3_files = sorted(Path(input_dir).glob("*.gff3"))
    
    if not gff3_files:
        print("No .gff3 files found!")
        return
    
    print(f"Found {len(gff3_files)} .gff3 files")
    
    # Step 2: Process all files
    all_results = []
    
    for gff3_file in tqdm(gff3_files, desc="Processing GFF3 files"):
        results = process_gff3_file(gff3_file)
        all_results.extend(results)
    
    # Step 3: Create DataFrame and save to CSV
    df = pd.DataFrame(all_results)
    df.to_csv(output_file, index=False)
    
    print(f"\nProcessing complete!")
    print(f"Total protein entries: {len(all_results)}")
    print(f"Output file: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Process GFF3 files to extract protein information. "
                    "Combines sequence ID (column 1) with protein ID from attributes (column 9) "
                    "and extracts start (column 4) and end (column 5) positions."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input directory containing .gff3 files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file path"
    )
    
    args = parser.parse_args()
    
    process_all_gff3_files(args.input, args.output)

if __name__ == "__main__":
    main()

