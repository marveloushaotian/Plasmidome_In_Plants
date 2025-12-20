#!/usr/bin/env python3
"""
Extract mapping between sample names and contig names from labeled FASTA files
"""

import argparse
import os
import sys
from pathlib import Path
import csv

def extract_sample_name(filename):
    """Extract sample name from filename (remove _labeled.fasta suffix)"""
    if filename.endswith('_labeled.fasta'):
        return filename[:-14]
    elif filename.endswith('.fasta'):
        return filename[:-6]
    else:
        return filename

def extract_contig_names(fasta_file):
    """Extract contig names from FASTA file (without '>' symbol)"""
    contig_names = []
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    contig_name = line[1:]  # Remove '>' symbol
                    contig_names.append(contig_name)
    except Exception as e:
        print(f"Warning: Error reading file {fasta_file}: {e}", file=sys.stderr)
    return contig_names

def process_labeled_directory(input_dir, output_file):
    """Process all labeled FASTA files in directory and create mapping table"""
    input_path = Path(input_dir)
    fasta_files = sorted(input_path.glob('*_labeled.fasta'))
    
    if not fasta_files:
        print(f"Warning: No *_labeled.fasta files found in {input_dir}", file=sys.stderr)
        return
    
    print(f"Found {len(fasta_files)} fasta files")
    
    all_mappings = []
    total = len(fasta_files)
    
    for idx, fasta_file in enumerate(fasta_files, 1):
        print(f"Processing file {idx}/{total}: {fasta_file.name}")
        
        sample_name = extract_sample_name(fasta_file.name)
        contig_names = extract_contig_names(fasta_file)
        
        for contig_name in contig_names:
            all_mappings.append({
                'sample_id': sample_name,
                'contig_name': contig_name
            })
    
    print(f"\nExtracted {len(all_mappings)} mappings")
    print(f"Writing to {output_file} ...")
    
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['sample_id', 'contig_name']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for mapping in all_mappings:
            writer.writerow(mapping)
    
    print(f"Done! Results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Extract mapping between sample names and contig names from labeled FASTA files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Extract mappings and output to CSV file
    python extract_contig_mapping.py -i Results/labeled -o contig_mapping.csv
    
Output format:
    CSV file with two columns:
    - sample_id: sample name (without _labeled.fasta suffix)
    - contig_name: contig name (without > symbol)
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input directory path (labeled directory)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output CSV file path'
    )
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.input):
        print(f"Error: Input directory does not exist: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    process_labeled_directory(args.input, args.output)

if __name__ == '__main__':
    main()

