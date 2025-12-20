#!/usr/bin/env python3
"""
Script to check paired-end FASTQ files for line count mismatches.

This script identifies paired FASTQ files (_1.fastq.gz and _2.fastq.gz)
and compares their line counts to find mismatches.

Usage:
    python check_paired_reads.py -i <input_dir> -o <output_file>
    python check_paired_reads.py -h  # Show help
"""

import os
import subprocess
import argparse
from pathlib import Path
from collections import defaultdict
import csv
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

def count_lines_in_gz(file_path):
    """
    Count lines in a gzipped file using zcat.
    
    Args:
        file_path: Path to the gzipped file
        
    Returns:
        tuple: (file_path, line_count) or (file_path, -1) if error
    """
    try:
        # Use zcat to decompress and wc -l to count lines
        result = subprocess.run(
            f"zcat {file_path} | wc -l",
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
        line_count = int(result.stdout.strip())
        return (file_path, line_count)
    except Exception as e:
        print(f"Error counting lines in {file_path}: {e}")
        return (file_path, -1)

def find_paired_files(input_dir):
    """
    Find all paired FASTQ files in the directory.
    
    Args:
        input_dir: Directory containing FASTQ files
        
    Returns:
        dict: Dictionary with sample names as keys and tuples of (file1, file2) as values
    """
    paired_files = defaultdict(lambda: [None, None])
    
    # Step 1: Find all fastq.gz files
    for file in Path(input_dir).glob("*.fastq.gz"):
        filename = file.name
        
        # Check if it's a _1 or _2 file
        if filename.endswith("_1.fastq.gz"):
            sample_name = filename[:-len("_1.fastq.gz")]
            paired_files[sample_name][0] = str(file)
        elif filename.endswith("_2.fastq.gz"):
            sample_name = filename[:-len("_2.fastq.gz")]
            paired_files[sample_name][1] = str(file)
    
    # Step 2: Filter out incomplete pairs
    complete_pairs = {}
    for sample_name, (file1, file2) in paired_files.items():
        if file1 and file2:
            complete_pairs[sample_name] = (file1, file2)
        else:
            print(f"Warning: Incomplete pair for sample {sample_name}")
    
    return complete_pairs

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Check paired-end FASTQ files for line count mismatches."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input directory containing paired FASTQ files"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file path"
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=8,
        help="Number of parallel jobs (default: 8)"
    )
    
    args = parser.parse_args()
    
    input_dir = args.input
    output_file = args.output
    num_jobs = args.jobs
    
    print(f"Scanning directory: {input_dir}")
    
    # Step 3: Find paired files
    paired_files = find_paired_files(input_dir)
    print(f"Found {len(paired_files)} paired samples")
    
    if len(paired_files) == 0:
        print("No paired files found. Exiting.")
        return
    
    # Step 4: Count lines in all files using parallel processing
    print("Counting lines in FASTQ files...")
    all_files = []
    for file1, file2 in paired_files.values():
        all_files.extend([file1, file2])
    
    line_counts = {}
    with ProcessPoolExecutor(max_workers=num_jobs) as executor:
        futures = {executor.submit(count_lines_in_gz, f): f for f in all_files}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
            file_path, count = future.result()
            line_counts[file_path] = count
    
    # Step 5: Compare paired files and identify mismatches
    print("\nComparing paired files...")
    results = []
    mismatches = []
    
    for sample_name, (file1, file2) in sorted(paired_files.items()):
        count1 = line_counts.get(file1, -1)
        count2 = line_counts.get(file2, -1)
        
        match_status = "Match" if count1 == count2 else "Mismatch"
        difference = abs(count1 - count2)
        
        result = {
            'Sample': sample_name,
            'File_R1': os.path.basename(file1),
            'File_R2': os.path.basename(file2),
            'Lines_R1': count1,
            'Lines_R2': count2,
            'Difference': difference,
            'Status': match_status
        }
        
        results.append(result)
        if match_status == "Mismatch":
            mismatches.append(result)
    
    # Step 6: Write results to CSV
    print(f"\nWriting results to {output_file}...")
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Sample', 'File_R1', 'File_R2', 'Lines_R1', 'Lines_R2', 'Difference', 'Status']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for result in results:
            writer.writerow(result)
    
    # Step 7: Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total pairs checked: {len(results)}")
    print(f"Matching pairs: {len(results) - len(mismatches)}")
    print(f"Mismatched pairs: {len(mismatches)}")
    
    if mismatches:
        print("\nMismatched samples:")
        for m in mismatches:
            print(f"  - {m['Sample']}: R1={m['Lines_R1']}, R2={m['Lines_R2']}, Diff={m['Difference']}")
    
    print(f"\nResults saved to: {output_file}")

if __name__ == "__main__":
    main()

