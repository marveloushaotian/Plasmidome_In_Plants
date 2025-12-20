#!/usr/bin/env python3
"""
Statistical analysis of genus-level defense systems, anti-defense systems, and AMR.

This script processes contig-level data to generate genus-level statistics for:
1. Genus counts across different contig types
2. Top 20 defense subtypes
3. All anti-defense system types
4. Top 20 AMR types

Usage:
    python 205_genus_statistics.py -g <genus_file> -i <input_file> -o <output_prefix> [-h]

Arguments:
    -g: Input CSV file containing genus names (GTDB_Genus column)
    -i: Input CSV file with contig-level annotations
    -o: Output file prefix for results
    -h: Show this help message

Example:
    python 205_genus_statistics.py -g Result/NCBI_4395_Batch/05_Tree/Genus_Level/Genus_Name.csv -i Result/NCBI_4395_Batch/Contig_Sample_Mapping_Final_with_Provirus_Overlap_GTDB_corrected.csv -o genus_stats
"""

import pandas as pd
import argparse
import numpy as np
from collections import defaultdict, Counter
from tqdm import tqdm


def parse_defense_subtypes(subtype_str):
    """
    Parse Defense_Subtype string and count occurrences.
    
    Args:
        subtype_str: String containing defense subtypes (may be comma-separated)
        
    Returns:
        List of individual defense subtypes
    """
    if pd.isna(subtype_str) or subtype_str == '':
        return []
    
    # Split by comma and strip whitespace
    subtypes = [s.strip() for s in str(subtype_str).split(',') if s.strip()]
    return subtypes


def should_exclude_defense(defense_name):
    """
    Check if a defense should be excluded based on filtering criteria.
    
    Args:
        defense_name: Defense system name
        
    Returns:
        True if should be excluded, False otherwise
    """
    # Exclude specific defenses
    exclude_list = ["VSPR", "dXTPase", "HEC-01", "HEC-09", "PifA"]
    if defense_name in exclude_list:
        return True
    
    # Exclude defenses containing "other"
    if 'other' in defense_name.lower():
        return True
    
    return False


def determine_contig_type(row):
    """
    Determine the final contig type considering Provirus_Overlap.
    
    Args:
        row: DataFrame row
        
    Returns:
        String indicating contig type (Chromosome, Plasmid, or Virus)
    """
    contig_type2 = row['Contig_Type2']
    provirus_overlap = row['Provirus_Overlap']
    
    # Default type
    final_type = contig_type2
    
    # Check for Defense_on_Provirus or AMR_on_Provirus
    if pd.notna(provirus_overlap) and provirus_overlap != '':
        overlap_items = [s.strip() for s in str(provirus_overlap).split(',')]
        
        # If Defense_on_Provirus or AMR_on_Provirus exists and original type is Chromosome
        if contig_type2 == 'Chromosome':
            if 'Defense_on_Provirus' in overlap_items or 'AMR_on_Provirus' in overlap_items:
                final_type = 'Virus'
    
    return final_type


def complete_genus_list(df, genus_list, genus_column='Genus', value_columns=None):
    """
    Ensure DataFrame contains all genera from genus_list, filling missing ones with 0.
    
    Args:
        df: Input DataFrame
        genus_list: Complete list of genus names in desired order
        genus_column: Name of the genus column
        value_columns: List of value columns to fill with 0 for missing genera
        
    Returns:
        Complete DataFrame with all genera
    """
    # Create a complete DataFrame with all genera
    complete_df = pd.DataFrame({genus_column: genus_list})
    
    # Merge with existing data
    merged_df = pd.merge(complete_df, df, on=genus_column, how='left')
    
    # Fill NaN values with 0 for value columns
    if value_columns:
        for col in value_columns:
            if col in merged_df.columns:
                merged_df[col] = merged_df[col].fillna(0).astype(int)
    else:
        # Fill all numeric columns with 0
        numeric_cols = merged_df.select_dtypes(include=['number']).columns
        merged_df[numeric_cols] = merged_df[numeric_cols].fillna(0).astype(int)
    
    return merged_df


def calculate_mean_plasmid_percent_by_sample(df, genus_list):
    """
    Calculate mean plasmid_percent per genus based on unique samples.
    
    For each Host-Genus combination, take one plasmid_percent value per Sample_ID,
    then calculate the mean across all samples.
    
    Args:
        df: Input DataFrame
        genus_list: List of genus names to analyze (in order)
        
    Returns:
        Dictionary with Host as keys and DataFrames with Genus and Mean_Plasmid_Percent
    """
    print("Calculating mean plasmid percent by sample...")
    
    # Filter for target genera
    df_filtered = df[df['Genus_CRBC_Updated'].isin(genus_list)].copy()
    
    # For each Host and Genus, get unique Sample_ID and their plasmid_percent
    result_dict = {}
    
    for host in df_filtered['Host'].unique():
        host_data = df_filtered[df_filtered['Host'] == host]
        
        # For each genus, calculate mean plasmid_percent across unique samples
        genus_plasmid = []
        
        for genus in genus_list:
            genus_data = host_data[host_data['Genus_CRBC_Updated'] == genus]
            
            if len(genus_data) > 0:
                # Get one plasmid_percent per Sample_ID (take first occurrence)
                sample_plasmid = genus_data.groupby('Sample_ID')['plasmid_percent'].first()
                mean_plasmid = sample_plasmid.mean()
            else:
                mean_plasmid = 0.0
            
            genus_plasmid.append({
                'Genus': genus,
                'Mean_Plasmid_Percent': mean_plasmid
            })
        
        result_dict[host] = pd.DataFrame(genus_plasmid)
    
    return result_dict


def count_genus_numbers_and_plasmid_percent(df, genus_list):
    """
    Count number of unique Sample_IDs per genus (by Host only) and calculate mean plasmid_percent.
    Also calculate total length (in kb) per genus/host/contig_type based on unique samples.
    
    Note: Genus_Number is now calculated based on unique Sample_IDs per Host, not per Contig_Type.
    This means each Host has a fixed Genus_Number regardless of Contig_Type.
    
    Args:
        df: Input DataFrame
        genus_list: List of genus names to analyze (in order)
        
    Returns:
        Dictionary with (Host, Contig_Type) as keys and DataFrames as values
    """
    print("Step 1: Counting genus numbers (unique Sample_IDs per Host) and calculating total lengths...")
    
    # Add final contig type column
    df['Final_Contig_Type'] = df.apply(determine_contig_type, axis=1)
    
    # Filter for target genera
    df_filtered = df[df['Genus_CRBC_Updated'].isin(genus_list)].copy()
    
    # Create a global taxonomy dictionary for all genera in the dataset
    taxonomy_cols = ['Kingdom_CRBC', 'Phylum_CRBC', 'Class_CRBC', 'Order_CRBC', 'Family_CRBC']
    global_taxonomy = df_filtered.groupby('Genus_CRBC_Updated')[taxonomy_cols].first().reset_index()
    global_taxonomy.columns = ['Genus', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']
    global_taxonomy_dict = global_taxonomy.set_index('Genus').to_dict('index')
    
    # Calculate mean plasmid percent by sample (per Host only, not by Contig_Type)
    plasmid_percent_dict = calculate_mean_plasmid_percent_by_sample(df, genus_list)
    
    # NEW: Calculate Genus_Number based on unique Sample_IDs per Host
    # Also extract taxonomy information for each genus
    # Group by Host and Genus, then count unique Sample_IDs
    genus_number_dict = {}
    genus_taxonomy_dict = {}
    
    for host, host_data in df_filtered.groupby('Host'):
        # Count unique Sample_IDs per genus for this host
        genus_counts = host_data.groupby('Genus_CRBC_Updated')['Sample_ID'].nunique().reset_index()
        genus_counts.columns = ['Genus', 'Genus_Number']
        genus_number_dict[host] = genus_counts
        
        # Get taxonomy information for each genus (take first occurrence)
        genus_taxonomy = host_data.groupby('Genus_CRBC_Updated')[taxonomy_cols].first().reset_index()
        genus_taxonomy.columns = ['Genus', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']
        genus_taxonomy_dict[host] = genus_taxonomy
    
    # NEW: Calculate total length per genus/host/contig_type
    # For each genus/host/contig_type, sum the length from unique samples
    print("  Calculating total lengths per genus/host/contig_type...")
    length_dict = {}
    
    for (host, contig_type), group in df_filtered.groupby(['Host', 'Final_Contig_Type']):
        genus_lengths = []
        
        for genus in genus_list:
            genus_data = group[group['Genus_CRBC_Updated'] == genus]
            
            if len(genus_data) > 0:
                # Get unique samples for this genus
                unique_samples = genus_data['Sample_ID'].unique()
                
                # For each unique sample, get the corresponding length once
                total_length = 0
                for sample_id in unique_samples:
                    sample_row = genus_data[genus_data['Sample_ID'] == sample_id].iloc[0]
                    
                    # Select the appropriate length based on contig type
                    if contig_type == 'Chromosome':
                        total_length += sample_row['chromosome_length']
                    elif contig_type == 'Plasmid':
                        total_length += sample_row['plasmid_length']
                    elif contig_type == 'Virus':
                        total_length += sample_row['virus_length']
                
                # Convert to kb
                total_length_kb = total_length / 1000.0
            else:
                total_length_kb = 0.0
            
            genus_lengths.append({
                'Genus': genus,
                'Total_Length_kb': total_length_kb
            })
        
        length_dict[(host, contig_type)] = pd.DataFrame(genus_lengths)
    
    # Get all unique Host and Final_Contig_Type combinations
    contig_type_combinations = df_filtered[['Host', 'Final_Contig_Type']].drop_duplicates()
    
    # Split by Host and Contig_Type
    result_dict = {}
    for _, row in contig_type_combinations.iterrows():
        host = row['Host']
        contig_type = row['Final_Contig_Type']
        
        # Get Genus_Number for this host (same for all contig types)
        if host in genus_number_dict:
            group_data = genus_number_dict[host].copy()
        else:
            group_data = pd.DataFrame({'Genus': [], 'Genus_Number': []})
        
        # Merge with taxonomy data for this host
        if host in genus_taxonomy_dict:
            group_data = pd.merge(group_data, genus_taxonomy_dict[host], on='Genus', how='left')
        
        # Merge with plasmid percent data for this host
        if host in plasmid_percent_dict:
            group_data = pd.merge(group_data, plasmid_percent_dict[host], on='Genus', how='left')
            group_data['Mean_Plasmid_Percent'] = group_data['Mean_Plasmid_Percent'].fillna(0.0)
        else:
            group_data['Mean_Plasmid_Percent'] = 0.0
        
        # Merge with length data for this host/contig_type combination
        if (host, contig_type) in length_dict:
            group_data = pd.merge(group_data, length_dict[(host, contig_type)], on='Genus', how='left')
            group_data['Total_Length_kb'] = group_data['Total_Length_kb'].fillna(0.0)
        else:
            group_data['Total_Length_kb'] = 0.0
        
        # Complete with all genera from genus_list
        # First complete the DataFrame
        complete_genus_df = pd.DataFrame({'Genus': genus_list})
        complete_data = pd.merge(complete_genus_df, group_data, on='Genus', how='left')
        
        # Fill NaN values: Genus_Number with 0, Mean_Plasmid_Percent with plasmid data if exists
        complete_data['Genus_Number'] = complete_data['Genus_Number'].fillna(0).astype(int)
        complete_data['Total_Length_kb'] = complete_data['Total_Length_kb'].fillna(0.0)
        
        # Fill taxonomy information from global_taxonomy_dict
        taxonomy_cols_to_fill = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family']
        for col in taxonomy_cols_to_fill:
            if col in complete_data.columns:
                for idx, row in complete_data.iterrows():
                    if pd.isna(row[col]) and row['Genus'] in global_taxonomy_dict:
                        complete_data.at[idx, col] = global_taxonomy_dict[row['Genus']][col]
        
        # For Mean_Plasmid_Percent, merge again with host plasmid data for missing genera
        if host in plasmid_percent_dict:
            # Update NaN values with data from plasmid_percent_dict
            for idx, row in complete_data.iterrows():
                if pd.isna(row['Mean_Plasmid_Percent']):
                    genus_plasmid = plasmid_percent_dict[host][plasmid_percent_dict[host]['Genus'] == row['Genus']]['Mean_Plasmid_Percent']
                    if len(genus_plasmid) > 0:
                        complete_data.at[idx, 'Mean_Plasmid_Percent'] = genus_plasmid.values[0]
                    else:
                        complete_data.at[idx, 'Mean_Plasmid_Percent'] = 0.0
        else:
            complete_data['Mean_Plasmid_Percent'] = complete_data['Mean_Plasmid_Percent'].fillna(0.0)
        
        result_dict[(host, contig_type)] = complete_data
    
    return result_dict


def count_defense_subtypes(df, genus_list, top_n=12):
    """
    Count defense subtypes per genus grouped by Host and final Contig_Type.
    Top N defenses are counted individually, rest are combined as 'Defense_Others'.
    
    Args:
        df: Input DataFrame with Final_Contig_Type column
        genus_list: List of genus names to analyze (in order)
        top_n: Number of top defenses to include
        
    Returns:
        Tuple of (top_defenses list, dictionary with results per group)
    """
    print("Step 2: Counting defense subtypes...")
    
    # Filter for target genera
    df_filtered = df[df['Genus_CRBC_Updated'].isin(genus_list)].copy()
    
    # Count all defenses to find top N
    all_defenses = []
    for _, row in tqdm(df_filtered.iterrows(), total=len(df_filtered), desc="Parsing defenses"):
        subtypes = parse_defense_subtypes(row['Defense_Subtype'])
        all_defenses.extend(subtypes)
    
    # Filter out excluded defenses
    filtered_defenses = [d for d in all_defenses if not should_exclude_defense(d)]
    
    # Get top N defenses
    defense_counter = Counter(filtered_defenses)
    top_defenses = [d[0] for d in defense_counter.most_common(top_n)]
    
    print(f"Top {top_n} defenses: {top_defenses}")
    
    # Count defenses per genus/host/type
    results = []
    
    grouped = df_filtered.groupby(['Host', 'Final_Contig_Type', 'Genus_CRBC_Updated'])
    
    for (host, contig_type, genus), group in tqdm(grouped, desc="Counting by genus"):
        defense_counts = Counter()
        others_count = 0
        
        for _, row in group.iterrows():
            subtypes = parse_defense_subtypes(row['Defense_Subtype'])
            for subtype in subtypes:
                if not should_exclude_defense(subtype):
                    if subtype in top_defenses:
                        defense_counts[subtype] += 1
                    else:
                        others_count += 1
        
        # Create result row
        result_row = {
            'Host': host,
            'Contig_Type': contig_type,
            'Genus': genus
        }
        for defense in top_defenses:
            result_row[defense] = defense_counts.get(defense, 0)
        result_row['Defense_Others'] = others_count
        
        results.append(result_row)
    
    df_results = pd.DataFrame(results)
    
    # Split by Host and Contig_Type
    result_dict = {}
    defense_cols_with_others = top_defenses + ['Defense_Others']
    for (host, contig_type), group in df_results.groupby(['Host', 'Contig_Type']):
        # Keep only Genus and defense columns
        cols_to_keep = ['Genus'] + defense_cols_with_others
        group_data = group[cols_to_keep]
        # Complete with all genera from genus_list
        complete_data = complete_genus_list(group_data, genus_list, 'Genus', defense_cols_with_others)
        result_dict[(host, contig_type)] = complete_data
    
    return top_defenses, result_dict


def count_antids_types(df, genus_list):
    """
    Count anti-defense system types per genus grouped by Host and final Contig_Type.
    
    Args:
        df: Input DataFrame with Final_Contig_Type column
        genus_list: List of genus names to analyze (in order)
        
    Returns:
        Tuple of (unique_antids list, dictionary with results per group)
    """
    print("Step 3: Counting anti-defense system types...")
    
    # Filter for target genera
    df_filtered = df[df['Genus_CRBC_Updated'].isin(genus_list)].copy()
    
    # Get all unique AntiDS types
    all_antids = []
    for val in df_filtered['AntiDS_Type'].dropna():
        if val != '':
            all_antids.extend([s.strip() for s in str(val).split(',')])
    
    unique_antids = sorted(list(set(all_antids)))
    print(f"Found {len(unique_antids)} unique anti-defense systems: {unique_antids}")
    
    # Count AntiDS per genus/host/type
    results = []
    
    grouped = df_filtered.groupby(['Host', 'Final_Contig_Type', 'Genus_CRBC_Updated'])
    
    for (host, contig_type, genus), group in tqdm(grouped, desc="Counting AntiDS by genus"):
        antids_counts = Counter()
        
        for _, row in group.iterrows():
            if pd.notna(row['AntiDS_Type']) and row['AntiDS_Type'] != '':
                types = [s.strip() for s in str(row['AntiDS_Type']).split(',')]
                for atype in types:
                    antids_counts[atype] += 1
        
        # Create result row
        result_row = {
            'Host': host,
            'Contig_Type': contig_type,
            'Genus': genus
        }
        for antids in unique_antids:
            result_row[antids] = antids_counts.get(antids, 0)
        
        results.append(result_row)
    
    df_results = pd.DataFrame(results)
    
    # Split by Host and Contig_Type
    result_dict = {}
    for (host, contig_type), group in df_results.groupby(['Host', 'Contig_Type']):
        # Keep only Genus and AntiDS columns
        cols_to_keep = ['Genus'] + unique_antids
        group_data = group[cols_to_keep]
        # Complete with all genera from genus_list
        complete_data = complete_genus_list(group_data, genus_list, 'Genus', unique_antids)
        result_dict[(host, contig_type)] = complete_data
    
    return unique_antids, result_dict


def count_amr_types(df, genus_list, top_n=12):
    """
    Count AMR types per genus grouped by Host and final Contig_Type.
    Top N AMR types are counted individually, rest are combined as 'AMR_Others'.
    Considers AMR_on_Provirus for Virus classification.
    
    Args:
        df: Input DataFrame with Final_Contig_Type column
        genus_list: List of genus names to analyze (in order)
        top_n: Number of top AMR types to include
        
    Returns:
        Tuple of (top_amr list, dictionary with results per group)
    """
    print("Step 4: Counting AMR types...")
    
    # Filter for target genera
    df_filtered = df[df['Genus_CRBC_Updated'].isin(genus_list)].copy()
    
    # For AMR, we need to recalculate Final_Contig_Type considering only AMR_on_Provirus
    def determine_amr_contig_type(row):
        contig_type2 = row['Contig_Type2']
        provirus_overlap = row['Provirus_Overlap']
        
        final_type = contig_type2
        
        if pd.notna(provirus_overlap) and provirus_overlap != '':
            overlap_items = [s.strip() for s in str(provirus_overlap).split(',')]
            if contig_type2 == 'Chromosome' and 'AMR_on_Provirus' in overlap_items:
                final_type = 'Virus'
        
        return final_type
    
    df_filtered['AMR_Contig_Type'] = df_filtered.apply(determine_amr_contig_type, axis=1)
    
    # Count all AMR types to find top N
    all_amr = []
    for val in df_filtered['AMR_Type'].dropna():
        if val != '':
            all_amr.extend([s.strip() for s in str(val).split(',')])
    
    amr_counter = Counter(all_amr)
    top_amr = [a[0] for a in amr_counter.most_common(top_n)]
    
    print(f"Top {top_n} AMR types: {top_amr}")
    
    # Count AMR per genus/host/type
    results = []
    
    grouped = df_filtered.groupby(['Host', 'AMR_Contig_Type', 'Genus_CRBC_Updated'])
    
    for (host, contig_type, genus), group in tqdm(grouped, desc="Counting AMR by genus"):
        amr_counts = Counter()
        others_count = 0
        
        for _, row in group.iterrows():
            if pd.notna(row['AMR_Type']) and row['AMR_Type'] != '':
                types = [s.strip() for s in str(row['AMR_Type']).split(',')]
                for atype in types:
                    if atype in top_amr:
                        amr_counts[atype] += 1
                    else:
                        others_count += 1
        
        # Create result row
        result_row = {
            'Host': host,
            'Contig_Type': contig_type,
            'Genus': genus
        }
        for amr in top_amr:
            result_row[amr] = amr_counts.get(amr, 0)
        result_row['AMR_Others'] = others_count
        
        results.append(result_row)
    
    df_results = pd.DataFrame(results)
    
    # Split by Host and Contig_Type
    result_dict = {}
    amr_cols_with_others = top_amr + ['AMR_Others']
    for (host, contig_type), group in df_results.groupby(['Host', 'Contig_Type']):
        # Keep only Genus and AMR columns
        cols_to_keep = ['Genus'] + amr_cols_with_others
        group_data = group[cols_to_keep]
        # Complete with all genera from genus_list
        complete_data = complete_genus_list(group_data, genus_list, 'Genus', amr_cols_with_others)
        result_dict[(host, contig_type)] = complete_data
    
    return top_amr, result_dict


def main():
    """Main function to parse arguments and execute analysis."""
    parser = argparse.ArgumentParser(
        description="Generate genus-level statistics for defense systems, anti-DS, and AMR",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '-g', '--genus',
        required=True,
        help='Input CSV file containing genus names (GTDB_Genus column)'
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input CSV file with contig-level annotations'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file prefix for results'
    )
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("Genus-Level Statistical Analysis")
    print("=" * 80)
    
    # Step 1: Load genus list
    print("\nLoading genus list...")
    genus_df = pd.read_csv(args.genus)
    
    # Get unique genus list from Genus_CRBC_Updated column
    genus_list = genus_df['Genus_CRBC_Updated'].unique().tolist()
    print(f"Loaded {len(genus_list)} unique genera")
    
    # Step 2: Load main data
    print("\nLoading contig data...")
    df = pd.read_csv(args.input, low_memory=False)
    print(f"Loaded {len(df)} rows")
    
    print("\n" + "=" * 80)
    
    # Step 3: Count genus numbers and calculate plasmid percent
    genus_counts_dict = count_genus_numbers_and_plasmid_percent(df, genus_list)
    print(f"Saving genus counts for {len(genus_counts_dict)} groups...")
    for (host, contig_type), data in genus_counts_dict.items():
        output_file = f"{args.output}_genus_counts_{host}_{contig_type}.csv"
        data.to_csv(output_file, index=False)
        print(f"  - Saved {output_file} ({len(data)} genera)")
    
    print("\n" + "=" * 80)
    
    # Step 4: Count defense subtypes (top 12 + Others)
    top_defenses, defense_counts_dict = count_defense_subtypes(df, genus_list, top_n=12)
    print(f"Saving defense counts for {len(defense_counts_dict)} groups...")
    for (host, contig_type), data in defense_counts_dict.items():
        output_file = f"{args.output}_defense_top12_{host}_{contig_type}.csv"
        data.to_csv(output_file, index=False)
        print(f"  - Saved {output_file} ({len(data)} genera)")
    
    print("\n" + "=" * 80)
    
    # Step 5: Count anti-defense systems (all)
    unique_antids, antids_counts_dict = count_antids_types(df, genus_list)
    print(f"Saving anti-DS counts for {len(antids_counts_dict)} groups...")
    for (host, contig_type), data in antids_counts_dict.items():
        output_file = f"{args.output}_antids_all_{host}_{contig_type}.csv"
        data.to_csv(output_file, index=False)
        print(f"  - Saved {output_file} ({len(data)} genera)")
    
    print("\n" + "=" * 80)
    
    # Step 6: Count AMR types (top 12 + Others)
    top_amr, amr_counts_dict = count_amr_types(df, genus_list, top_n=12)
    print(f"Saving AMR counts for {len(amr_counts_dict)} groups...")
    for (host, contig_type), data in amr_counts_dict.items():
        output_file = f"{args.output}_amr_top12_{host}_{contig_type}.csv"
        data.to_csv(output_file, index=False)
        print(f"  - Saved {output_file} ({len(data)} genera)")
    
    print("\n" + "=" * 80)
    print("Step 7: Merging all statistics and calculating derived metrics...")
    
    # Get all unique (Host, Contig_Type) combinations
    all_groups = set()
    all_groups.update(genus_counts_dict.keys())
    all_groups.update(defense_counts_dict.keys())
    all_groups.update(antids_counts_dict.keys())
    all_groups.update(amr_counts_dict.keys())
    
    print(f"Merging data for {len(all_groups)} groups...")
    
    for (host, contig_type) in sorted(all_groups):
        print(f"\n  Processing {host}_{contig_type}...")
        
        # Start with genus counts and plasmid percent
        if (host, contig_type) in genus_counts_dict:
            merged_df = genus_counts_dict[(host, contig_type)].copy()
            print(f"    - Added genus counts and plasmid percent ({len(merged_df)} genera)")
        else:
            print(f"    - Warning: No genus counts found")
            continue
        
        # Merge defense data
        defense_cols = []
        if (host, contig_type) in defense_counts_dict:
            defense_df = defense_counts_dict[(host, contig_type)]
            merged_df = pd.merge(merged_df, defense_df, on='Genus', how='left')
            # Fill NaN with 0 for defense counts
            defense_cols = [col for col in defense_df.columns if col != 'Genus']
            merged_df[defense_cols] = merged_df[defense_cols].fillna(0).astype(int)
            print(f"    - Merged defense systems ({len(defense_cols)} types)")
        
        # Merge anti-DS data
        antids_cols = []
        if (host, contig_type) in antids_counts_dict:
            antids_df = antids_counts_dict[(host, contig_type)]
            merged_df = pd.merge(merged_df, antids_df, on='Genus', how='left')
            # Fill NaN with 0 for anti-DS counts
            antids_cols = [col for col in antids_df.columns if col != 'Genus']
            merged_df[antids_cols] = merged_df[antids_cols].fillna(0).astype(int)
            print(f"    - Merged anti-defense systems ({len(antids_cols)} types)")
        
        # Merge AMR data
        amr_cols = []
        if (host, contig_type) in amr_counts_dict:
            amr_df = amr_counts_dict[(host, contig_type)]
            merged_df = pd.merge(merged_df, amr_df, on='Genus', how='left')
            # Fill NaN with 0 for AMR counts
            amr_cols = [col for col in amr_df.columns if col != 'Genus']
            merged_df[amr_cols] = merged_df[amr_cols].fillna(0).astype(int)
            print(f"    - Merged AMR types ({len(amr_cols)} types)")
        
        # Calculate derived metrics
        print(f"    - Calculating derived metrics...")
        
        # 1. Total defense count per genus and per kb length
        if defense_cols:
            # Include Defense_Others in defense_cols if it exists
            all_defense_cols = defense_cols.copy()
            if 'Defense_Others' in merged_df.columns and 'Defense_Others' not in all_defense_cols:
                all_defense_cols.append('Defense_Others')
            
            merged_df['Total_Defense_Count'] = merged_df[all_defense_cols].sum(axis=1)
            # Defense count per kb (normalized by total length in kb)
            merged_df['Defense_Per_kb'] = merged_df.apply(
                lambda row: row['Total_Defense_Count'] / row['Total_Length_kb'] if row['Total_Length_kb'] > 0 else 0,
                axis=1
            )
            # Individual defense per kb (including Defense_Others)
            for col in all_defense_cols:
                merged_df[f'{col}_Per_kb'] = merged_df.apply(
                    lambda row: row[col] / row['Total_Length_kb'] if row['Total_Length_kb'] > 0 else 0,
                    axis=1
                )
        
        # 2. Total anti-defense count per genus and per kb length
        if antids_cols:
            merged_df['Total_AntiDS_Count'] = merged_df[antids_cols].sum(axis=1)
            # Anti-defense count per kb (normalized by total length in kb)
            merged_df['AntiDS_Per_kb'] = merged_df.apply(
                lambda row: row['Total_AntiDS_Count'] / row['Total_Length_kb'] if row['Total_Length_kb'] > 0 else 0,
                axis=1
            )
            # Individual anti-defense per kb
            for col in antids_cols:
                merged_df[f'{col}_Per_kb'] = merged_df.apply(
                    lambda row: row[col] / row['Total_Length_kb'] if row['Total_Length_kb'] > 0 else 0,
                    axis=1
                )
        
        # 3. Total AMR count per genus and per kb length
        if amr_cols:
            # Include AMR_Others in amr_cols if it exists
            all_amr_cols = amr_cols.copy()
            if 'AMR_Others' in merged_df.columns and 'AMR_Others' not in all_amr_cols:
                all_amr_cols.append('AMR_Others')
            
            merged_df['Total_AMR_Count'] = merged_df[all_amr_cols].sum(axis=1)
            # AMR count per kb (normalized by total length in kb)
            merged_df['AMR_Per_kb'] = merged_df.apply(
                lambda row: row['Total_AMR_Count'] / row['Total_Length_kb'] if row['Total_Length_kb'] > 0 else 0,
                axis=1
            )
            # Individual AMR per kb (including AMR_Others)
            for col in all_amr_cols:
                merged_df[f'{col}_Per_kb'] = merged_df.apply(
                    lambda row: row[col] / row['Total_Length_kb'] if row['Total_Length_kb'] > 0 else 0,
                    axis=1
                )
        
        # Save merged file
        output_file = f"{args.output}_merged_{host}_{contig_type}.csv"
        merged_df.to_csv(output_file, index=False)
        print(f"    - Saved merged file: {output_file}")
        print(f"      ({len(merged_df)} genera × {len(merged_df.columns)} columns)")
    
    print("\n" + "=" * 80)
    print("Analysis completed successfully!")
    print(f"\nIndividual files generated:")
    print(f"  - Genus counts: {len(genus_counts_dict)} files")
    print(f"  - Defense systems: {len(defense_counts_dict)} files")
    print(f"  - Anti-defense systems: {len(antids_counts_dict)} files")
    print(f"  - AMR types: {len(amr_counts_dict)} files")
    print(f"\nMerged files generated: {len(all_groups)} files")
    print("=" * 80)


if __name__ == "__main__":
    main()

