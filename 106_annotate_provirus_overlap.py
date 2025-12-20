#!/usr/bin/env python3
"""
Annotate defense systems, PDC, AntiDS, and AMR that overlap with provirus regions

This script identifies provirus regions from the Contig_Type column 
(both Chromosome|Provirus and Plasmid|provirus) and checks if defense systems, 
PDC, AntiDS, or AMR regions overlap with any provirus regions.

Usage:
    python annotate_provirus_overlap.py -i input.csv -o output.csv
"""

import argparse
import pandas as pd
import re
from tqdm import tqdm


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Annotate features overlapping with provirus regions'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CSV file with contig-sample mapping'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output CSV file with provirus overlap annotations'
    )
    return parser.parse_args()


def extract_provirus_regions(contig_id, contig_type):
    """
    Extract provirus regions from Contig_ID column
    
    Args:
        contig_id: String like "NODE_3_length_625281_cov_88.590649|chromosome||provirus_region=513028-535250"
        contig_type: String like "Chromosome|Provirus", "Plasmid|Provirus", "Plasmid|provirus" to verify this is a provirus
    
    Returns:
        List of tuples [(start, end), ...] or empty list if no provirus regions
    """
    # Only process if contig_type indicates this is a provirus
    # Accept both "Chromosome|Provirus" and "Plasmid|Provirus" or "Plasmid|provirus" (case insensitive)
    contig_type_str = str(contig_type).lower()
    if pd.isna(contig_type) or 'provirus' not in contig_type_str:
        return []
    
    if pd.isna(contig_id) or 'provirus_region=' not in str(contig_id):
        return []
    
    regions = []
    # Find all provirus_region patterns
    pattern = r'provirus_region=(\d+)-(\d+)'
    matches = re.findall(pattern, str(contig_id))
    
    for match in matches:
        start = int(match[0])
        end = int(match[1])
        regions.append((start, end))
    
    return regions


def parse_positions(pos_string):
    """
    Parse position string that may contain multiple comma-separated values
    
    Args:
        pos_string: String like "1,650,491,758,615" or single value "513028"
    
    Returns:
        List of integers, or empty list if input is NaN or empty
    """
    if pd.isna(pos_string) or pos_string == '':
        return []
    
    # Convert to string and split by comma
    pos_str = str(pos_string).strip()
    if pos_str == '':
        return []
    
    try:
        positions = [int(x.strip()) for x in pos_str.split(',') if x.strip()]
        return positions
    except ValueError:
        return []


def check_overlap(start1, end1, start2, end2):
    """
    Check if two regions overlap
    
    Args:
        start1, end1: First region boundaries
        start2, end2: Second region boundaries
    
    Returns:
        True if regions overlap, False otherwise
    """
    # Two regions overlap if one starts before the other ends
    return not (end1 < start2 or end2 < start1)


def check_feature_overlap_with_provirus(start_positions, end_positions, provirus_regions):
    """
    Check if any feature region overlaps with any provirus region
    
    Args:
        start_positions: List of start positions
        end_positions: List of end positions
        provirus_regions: List of (start, end) tuples for provirus regions
    
    Returns:
        True if any feature region overlaps with any provirus region
    """
    if not start_positions or not end_positions or not provirus_regions:
        return False
    
    # Pair up start and end positions
    # If counts don't match, pair as many as possible
    num_features = min(len(start_positions), len(end_positions))
    
    for i in range(num_features):
        feature_start = start_positions[i]
        feature_end = end_positions[i]
        
        # Check against all provirus regions
        for prov_start, prov_end in provirus_regions:
            if check_overlap(feature_start, feature_end, prov_start, prov_end):
                return True
    
    return False


def annotate_provirus_overlaps(row, provirus_regions):
    """
    Annotate a row with provirus overlap information
    
    Args:
        row: DataFrame row
        provirus_regions: List of (start, end) tuples for provirus regions
    
    Returns:
        String with comma-separated annotations (e.g., "Defense_on_Provirus,AMR_on_Provirus")
    """
    annotations = []
    
    # Check Defense System overlap
    defense_starts = parse_positions(row['Defense_Start'])
    defense_ends = parse_positions(row['Defense_End'])
    if check_feature_overlap_with_provirus(defense_starts, defense_ends, provirus_regions):
        annotations.append('Defense_on_Provirus')
    
    # Check PDC overlap
    pdc_starts = parse_positions(row['PDC_Start'])
    pdc_ends = parse_positions(row['PDC_End'])
    if check_feature_overlap_with_provirus(pdc_starts, pdc_ends, provirus_regions):
        annotations.append('PDC_on_Provirus')
    
    # Check AntiDS overlap
    antids_starts = parse_positions(row['AntiDS_Start'])
    antids_ends = parse_positions(row['AntiDS_End'])
    if check_feature_overlap_with_provirus(antids_starts, antids_ends, provirus_regions):
        annotations.append('AntiDS_on_Provirus')
    
    # Check AMR overlap
    amr_starts = parse_positions(row['AMR_Start'])
    amr_ends = parse_positions(row['AMR_End'])
    if check_feature_overlap_with_provirus(amr_starts, amr_ends, provirus_regions):
        annotations.append('AMR_on_Provirus')
    
    return ','.join(annotations) if annotations else ''


def main():
    """Main function"""
    args = parse_args()
    
    # Step 1: Read input CSV
    print(f"Reading input file: {args.input}")
    df = pd.read_csv(args.input)
    print(f"Total rows: {len(df)}")
    
    # Step 2: Extract provirus regions for each row
    print("Extracting provirus regions...")
    provirus_regions_list = []
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing rows"):
        regions = extract_provirus_regions(row['Contig_ID'], row['Contig_Type'])
        provirus_regions_list.append(regions)
    
    # Step 3: Check for overlaps and annotate
    print("Checking for overlaps with provirus regions...")
    annotations = []
    for idx, (_, row) in enumerate(tqdm(df.iterrows(), total=len(df), desc="Annotating")):
        annotation = annotate_provirus_overlaps(row, provirus_regions_list[idx])
        annotations.append(annotation)
    
    # Step 4: Add annotation column to dataframe
    df['Provirus_Overlap'] = annotations
    
    # Step 5: Save output
    print(f"Saving results to: {args.output}")
    df.to_csv(args.output, index=False)
    
    # Step 6: Print summary statistics
    print("\n=== Summary ===")
    total_with_provirus = sum(1 for regions in provirus_regions_list if regions)
    total_with_overlap = sum(1 for ann in annotations if ann)
    
    print(f"Rows with provirus regions: {total_with_provirus}")
    print(f"Rows with overlapping features: {total_with_overlap}")
    
    # Count each type of overlap
    if total_with_overlap > 0:
        defense_count = sum(1 for ann in annotations if 'Defense_on_Provirus' in ann)
        pdc_count = sum(1 for ann in annotations if 'PDC_on_Provirus' in ann)
        antids_count = sum(1 for ann in annotations if 'AntiDS_on_Provirus' in ann)
        amr_count = sum(1 for ann in annotations if 'AMR_on_Provirus' in ann)
        
        print(f"\nOverlap counts:")
        print(f"  Defense_on_Provirus: {defense_count}")
        print(f"  PDC_on_Provirus: {pdc_count}")
        print(f"  AntiDS_on_Provirus: {antids_count}")
        print(f"  AMR_on_Provirus: {amr_count}")
    
    print("\nDone!")


if __name__ == '__main__':
    main()

