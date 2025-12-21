#!/usr/bin/env python3
"""
Merge defense system annotations from PADLOC, DefenseFinder, and CCTyper
COMPLETE OPTIMIZED VERSION - with all missing data cleaning steps
"""

import pandas as pd
import sys
from collections import defaultdict
import ast
import numpy as np
from tqdm import tqdm

# Enable tqdm for pandas operations
tqdm.pandas()

def load_mapping_file(mapping_file):
    """Load the defense systems name mapping file"""
    try:
        return pd.read_excel(mapping_file)
    except Exception as e:
        print(f"Error loading mapping file: {e}", file=sys.stderr)
        sys.exit(1)

def create_mapping_dicts(defense_systems_mapping):
    """Create mapping dictionaries for each tool"""
    mapping_dicts = {
        'defense_finder_type': defense_systems_mapping.set_index('deffinder_type')['unified type'].to_dict(),
        'defense_finder_subtype': defense_systems_mapping.set_index('deffinder_subtype')['unified subtype'].to_dict(),
        'padloc_type': defense_systems_mapping.set_index('padloc_subtype')['unified type'].to_dict(),
        'padloc_subtype': defense_systems_mapping.set_index('padloc_subtype')['unified subtype'].to_dict(),
        'cctyper_type': defense_systems_mapping.set_index('cctyper_subtype')['unified type'].to_dict(),
        'cctyper_subtype': defense_systems_mapping.set_index('cctyper_subtype')['unified subtype'].to_dict(),
    }
    return mapping_dicts

def process_defense_finder(df, mapping_dicts):
    """Process DefenseFinder data - OPTIMIZED"""
    print("[INFO] Processing DefenseFinder: applying unified naming...", file=sys.stderr)
    # Apply unified naming
    df['unified type'] = df['type'].replace(mapping_dicts['defense_finder_type'])
    df['unified subtype'] = df['subtype'].replace(mapping_dicts['defense_finder_subtype'])
    df['sys_id'] = df['sys_id'] + "|" + df["sys_beg"].astype(str).str.rsplit('_', n=1).str[0]

    print("[INFO] Processing DefenseFinder: splitting protein lists...", file=sys.stderr)
    # Split protein lists
    df['protein_in_syst'] = df['protein_in_syst'].str.split(',')
    df['name_of_profiles_in_sys'] = df['name_of_profiles_in_sys'].str.split(',')

    print("[INFO] Processing DefenseFinder: filtering multi-contig systems...", file=sys.stderr)
    # Filter multi-contig systems - VECTORIZED
    df['contigs_in_sys'] = df['protein_in_syst'].apply(
        lambda x: len(set([item.rsplit('_',1)[0] for item in x]))
    )
    df = df[df['contigs_in_sys'] == 1].copy()
    df.drop('contigs_in_sys', axis=1, inplace=True)

    print(f"[INFO] Processing DefenseFinder: expanding {len(df)} systems to protein level...", file=sys.stderr)
    # Expand to one row per protein - OPTIMIZED
    expanded_rows = []
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Expanding proteins"):
        for protein, profile in zip(row['protein_in_syst'], row['name_of_profiles_in_sys']):
            new_row = row.copy()
            new_row['protein_in_syst'] = protein
            new_row['name_of_profiles_in_sys'] = profile
            expanded_rows.append(new_row)
    
    expanded_df = pd.DataFrame(expanded_rows)

    print("[INFO] Processing DefenseFinder: splitting fragmented systems...", file=sys.stderr)
    # Split fragmented systems based on protein ordering - OPTIMIZED
    expanded_df['seqid'] = expanded_df['sys_id'].str.split('|').str[1]
    expanded_df['prot_ord'] = expanded_df['protein_in_syst'].astype(str).str.rsplit('_',n=1).str[1].astype(int)
    expanded_df = expanded_df.sort_values(['sys_id','prot_ord']).reset_index(drop=True)

    # Vectorized gap detection
    expanded_df['sys_id_shifted'] = expanded_df['sys_id'].shift(1)
    expanded_df['prot_ord_shifted'] = expanded_df['prot_ord'].shift(1)
    expanded_df['gap'] = expanded_df['prot_ord'] - expanded_df['prot_ord_shifted']
    expanded_df['is_new_sys'] = (expanded_df['sys_id'] != expanded_df['sys_id_shifted']) | (expanded_df['gap'] > 5)
    expanded_df['sub_sys_num'] = expanded_df.groupby('sys_id')['is_new_sys'].cumsum() - 1
    
    # Update sys_id with sub-system number where needed
    mask = expanded_df['sub_sys_num'] > 0
    expanded_df.loc[mask, 'sys_id'] = expanded_df.loc[mask, 'sys_id'] + "|" + expanded_df.loc[mask, 'sub_sys_num'].astype(str)

    expanded_df.drop(['seqid','prot_ord', 'sys_id_shifted', 'prot_ord_shifted', 'gap', 'is_new_sys', 'sub_sys_num'], axis=1, inplace=True)
    
    print(f"[INFO] DefenseFinder processing complete: {len(expanded_df)} protein entries", file=sys.stderr)
    return expanded_df

def process_padloc(df, mapping_dicts):
    """Process PADLOC data"""
    print("[INFO] Processing PADLOC data...", file=sys.stderr)
    df['unified type'] = df['system'].replace(mapping_dicts['padloc_type'])
    df['unified subtype'] = df['system'].replace(mapping_dicts['padloc_subtype'])
    df["system.number"] = df['seqid'] + "|" + df["system.number"].astype(str)
    print(f"[INFO] PADLOC processing complete: {len(df)} entries", file=sys.stderr)
    return df

def add_prot_ids_to_cctyper(df):
    """
    Add Prot_IDs column to CCTyper dataframe if it doesn't exist - OPTIMIZED
    """
    if 'Prot_IDs' in df.columns:
        print("[INFO] Prot_IDs column already exists in CCTyper data", file=sys.stderr)
        return df
    
    print("[INFO] Adding Prot_IDs column to CCTyper data...", file=sys.stderr)
    
    def create_prot_id(row):
        contig = row['Contig']
        if isinstance(row['Positions'], str):
            try:
                positions = ast.literal_eval(row['Positions'])
            except:
                positions = []
        else:
            positions = row['Positions']
        
        if isinstance(positions, list) and len(positions) > 0:
            prot_ids = [f"{contig}_{pos}" for pos in positions]
            return str(prot_ids)
        else:
            return "[]"
    
    # Vectorized operation with progress bar
    df['Prot_IDs'] = df.progress_apply(create_prot_id, axis=1)
    
    return df

def process_cctyper(df, mapping_dicts):
    """Process CCTyper data - OPTIMIZED"""
    print("[INFO] Processing CCTyper data...", file=sys.stderr)
    
    # Add Prot_IDs column if it doesn't exist
    df = add_prot_ids_to_cctyper(df)
    
    # Use Prediction column for mapping
    df['unified type'] = df['Prediction'].replace(mapping_dicts['cctyper_type'])
    df['unified subtype'] = df['Prediction'].replace(mapping_dicts['cctyper_subtype'])
    
    print("[INFO] CCTyper: converting lists...", file=sys.stderr)
    # Convert string representations to lists
    df['Genes'] = df['Genes'].apply(ast.literal_eval)
    df['Prot_IDs'] = df['Prot_IDs'].apply(ast.literal_eval)

    print(f"[INFO] CCTyper: expanding {len(df)} operons to protein level...", file=sys.stderr)
    # Expand to one row per protein - OPTIMIZED
    expanded_rows = []
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Expanding CCTyper"):
        for gene, prot_id in zip(row['Genes'], row['Prot_IDs']):
            new_row = row.copy()
            new_row['Genes'] = gene
            new_row['Prot_IDs'] = prot_id
            expanded_rows.append(new_row)
    
    expanded_df = pd.DataFrame(expanded_rows)
    
    print(f"[INFO] CCTyper processing complete: {len(expanded_df)} protein entries", file=sys.stderr)
    return expanded_df

def aggregate_merge(main_df, merge_df, left_on, right_on, new_columns):
    """Merge with aggregation - OPTIMIZED"""
    aggregated_dataframes = []

    for col in new_columns:
        aggregated_col = merge_df.groupby(right_on)[col].apply(
            lambda x: ','.join(x.dropna().unique())
        ).reset_index()
        aggregated_dataframes.append(aggregated_col)

    aggregated_df = pd.concat(aggregated_dataframes, axis=1)
    aggregated_df = aggregated_df.loc[:,~aggregated_df.columns.duplicated()]

    return main_df.merge(aggregated_df, left_on=left_on, right_on=right_on, how='left')

def select_preferred_value(values, preference_order):
    """Select preferred value based on tool priority"""
    for preference in preference_order:
        if preference in values:
            return preference
    return None if not values else values[0]

def determine_merged_info(prot_id_group, prot_id_df):
    """Determine merged type, subtype, and operon for a group of proteins"""
    # Get the info for all proteins in the group
    group_info = prot_id_df.loc[list(prot_id_group)]
    
    # Collect unique types, subtypes, and operons from each source
    unique_types = []
    unique_subtypes = []
    operon_names = {'cctyper': set(), 'deffinder': set(), 'padloc': set()}
    
    for prot_id, row in group_info.iterrows():
        for source in ['cctyper', 'deffinder', 'padloc']:
            # Handle comma-separated values (from groupby aggregation)
            if pd.notnull(row[source]['type']):
                for z in str(row[source]['type']).split(','):
                    if z and z not in unique_types:
                        unique_types.append(z)
            if pd.notnull(row[source]['subtype']):
                for z in str(row[source]['subtype']).split(','):
                    if z and z not in unique_subtypes:
                        unique_subtypes.append(z)
            if pd.notnull(row[source]['operon']):
                operon_names[source].add(str(row[source]['operon']))
    
    # Select the merged operon name based on the preference order
    merged_operon = None
    preference_order = ['cctyper', 'deffinder', 'padloc']
    for preference in preference_order:
        if operon_names[preference]:
            merged_operon = list(operon_names[preference])[0]
            break
    
    # Remove DMS and DRT from the unique types and subtypes if there are other options available
    if len(unique_types) > 1 and ('DMS' in unique_types or 'DRT' in unique_types):
        if 'DMS' in unique_types:
            unique_types.remove('DMS')
        if 'DRT' in unique_types:
            unique_types.remove('DRT')
    if len(unique_subtypes) > 1 and ('DMS' in unique_subtypes or 'DRT' in unique_subtypes):
        if 'DMS' in unique_subtypes:
            unique_subtypes.remove('DMS')
        if 'DRT' in unique_subtypes:
            unique_subtypes.remove('DRT')
    
    # Select the merged type and subtype based on the preference order of the sources
    merged_type = select_preferred_value(list(unique_types), preference_order)
    merged_subtype = select_preferred_value(list(unique_subtypes), preference_order)
    
    return pd.Series({'merged_type': merged_type, 'merged_subtype': merged_subtype, 'merged_operon': merged_operon})

def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Merge defense system annotations (COMPLETE OPTIMIZED VERSION)')
    parser.add_argument('--mapping', required=True, help='Defense systems name mapping file')
    parser.add_argument('--gff', required=True, help='GFF annotation file')
    parser.add_argument('--padloc', required=True, help='PADLOC results')
    parser.add_argument('--defensefinder', required=True, help='DefenseFinder results')
    parser.add_argument('--cctyper', required=True, help='CCTyper results')
    
    args = parser.parse_args()
    
    print("="*80, file=sys.stderr)
    print("[INFO] Starting merge process (COMPLETE OPTIMIZED VERSION)", file=sys.stderr)
    print("="*80, file=sys.stderr)
    
    print("[INFO] Loading input files...", file=sys.stderr)
    
    # Load mapping file
    defense_systems_mapping = load_mapping_file(args.mapping)
    mapping_dicts = create_mapping_dicts(defense_systems_mapping)
    
    # Load GFF file
    print("[INFO] Loading GFF file (this may take a while for large files)...", file=sys.stderr)
    all_gff_df = pd.read_csv(
        args.gff, sep='\t', header=0,
        names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "prot_id"]
    )
    print(f"[INFO] Loaded {len(all_gff_df):,} rows from GFF", file=sys.stderr)
    
    # Load tool results
    print("[INFO] Loading PADLOC results...", file=sys.stderr)
    padloc_df = pd.read_csv(args.padloc)  # CSV format
    print(f"[INFO] Loaded {len(padloc_df):,} PADLOC entries", file=sys.stderr)
    
    print("[INFO] Loading DefenseFinder results...", file=sys.stderr)
    defense_finder_df = pd.read_csv(args.defensefinder, sep='\t')
    print(f"[INFO] Loaded {len(defense_finder_df):,} DefenseFinder systems", file=sys.stderr)
    
    print("[INFO] Loading CCTyper results...", file=sys.stderr)
    cas_operons_df = pd.read_csv(args.cctyper, sep="\t")
    print(f"[INFO] Loaded {len(cas_operons_df):,} CCTyper operons", file=sys.stderr)
    
    print("\n" + "="*80, file=sys.stderr)
    print("[PROCESSING STAGE 1] Processing tool-specific data", file=sys.stderr)
    print("="*80, file=sys.stderr)
    
    expanded_defense_finder_df = process_defense_finder(defense_finder_df, mapping_dicts)
    padloc_df = process_padloc(padloc_df, mapping_dicts)
    expanded_cas_operons_df = process_cctyper(cas_operons_df, mapping_dicts)
    
    print("\n" + "="*80, file=sys.stderr)
    print("[PROCESSING STAGE 2] Merging all annotations", file=sys.stderr)
    print("="*80, file=sys.stderr)
    
    merged_gff_df = all_gff_df.copy()
    
    # Merge with expanded DefenseFinder data frame
    print("[INFO] Merging DefenseFinder annotations...", file=sys.stderr)
    merged_gff_df = aggregate_merge(
        merged_gff_df, 
        expanded_defense_finder_df, 
        'prot_id', 
        'protein_in_syst', 
        ['sys_id', 'unified type', 'unified subtype']
    )
    merged_gff_df.rename(columns={
        'sys_id': 'operon deffinder',
        'unified type': 'type deffinder',
        'unified subtype': 'subtype deffinder'
    }, inplace=True)
    
    # Merge with PADLOC data frame
    print("[INFO] Merging PADLOC annotations...", file=sys.stderr)
    merged_gff_df = aggregate_merge(
        merged_gff_df, 
        padloc_df, 
        'prot_id', 
        'target.name', 
        ['system.number', 'unified type', 'unified subtype']
    )
    merged_gff_df.rename(columns={
        'system.number': 'operon padloc',
        'unified type': 'type padloc',
        'unified subtype': 'subtype padloc'
    }, inplace=True)
    
    # Merge with expanded CAS operons data frame
    print("[INFO] Merging CCTyper annotations...", file=sys.stderr)
    merged_gff_df = aggregate_merge(
        merged_gff_df, 
        expanded_cas_operons_df, 
        'prot_id', 
        'Prot_IDs', 
        ['Operon', 'unified type', 'unified subtype']
    )
    merged_gff_df.rename(columns={
        'Operon': 'operon cctyper',
        'unified type': 'type cctyper',
        'unified subtype': 'subtype cctyper'
    }, inplace=True)
    
    # ===== CRITICAL FIX 1: Add groupby aggregation step =====
    print("[INFO] Aggregating duplicate protein IDs (CRITICAL STEP)...", file=sys.stderr)
    print(f"[INFO] Before aggregation: {len(merged_gff_df)} rows", file=sys.stderr)
    
    # Drop rows where all annotation columns are NaN
    merged_gff_df = merged_gff_df.dropna(
        subset=['type deffinder', 'subtype deffinder', 'operon deffinder',
                'type padloc', 'subtype padloc', 'operon padloc',
                'type cctyper', 'subtype cctyper', 'operon cctyper'],
        how='all'
    )
    
    # Group by prot_id and aggregate
    merged_gff_df = merged_gff_df.groupby('prot_id', as_index=False).agg({
        'seqid': lambda x: ','.join(map(str, x.dropna().unique())),
        'source': lambda x: ','.join(map(str, x.dropna().unique())),
        'type': lambda x: ','.join(map(str, x.dropna().unique())),
        'start': 'min',
        'end': 'max',
        'score': lambda x: ','.join(map(str, x.dropna().unique())),
        'strand': lambda x: ','.join(map(str, x.dropna().unique())),
        'phase': lambda x: ','.join(map(str, x.dropna().unique())) if 'phase' in x.index else '',
        'type deffinder': lambda x: ','.join(map(str, x.dropna().unique())),
        'subtype deffinder': lambda x: ','.join(map(str, x.dropna().unique())),
        'operon deffinder': lambda x: ','.join(map(str, x.dropna().unique())),
        'type padloc': lambda x: ','.join(map(str, x.dropna().unique())),
        'subtype padloc': lambda x: ','.join(map(str, x.dropna().unique())),
        'operon padloc': lambda x: ','.join(map(str, x.dropna().unique())),
        'type cctyper': lambda x: ','.join(map(str, x.dropna().unique())),
        'subtype cctyper': lambda x: ','.join(map(str, x.dropna().unique())),
        'operon cctyper': lambda x: ','.join(map(str, x.dropna().unique())),
    })
    
    print(f"[INFO] After aggregation: {len(merged_gff_df)} rows", file=sys.stderr)
    
    # ===== CRITICAL FIX 2: Split comma-separated Padloc values =====
    print("[INFO] Cleaning Padloc annotations (taking first value if multiple)...", file=sys.stderr)
    merged_gff_df['type padloc'] = merged_gff_df['type padloc'].str.split(',').str[0]
    merged_gff_df['subtype padloc'] = merged_gff_df['subtype padloc'].str.split(',').str[0]
    
    # Save intermediate merged file
    print("[INFO] Saving intermediate merged file (all.merged)...", file=sys.stderr)
    merged_gff_df.to_csv('all.merged', sep='\t', index=False)
    
    print("\n" + "="*80, file=sys.stderr)
    print("[PROCESSING STAGE 3] Resolving conflicts and creating unified annotations", file=sys.stderr)
    print("="*80, file=sys.stderr)
    
    # Create protein ID dictionary - OPTIMIZED using vectorized operations
    print("[INFO] Building protein ID dictionary...", file=sys.stderr)
    
    # Create a structured dictionary more efficiently
    prot_id_dict = {}
    for prot_id in tqdm(merged_gff_df['prot_id'].unique(), desc="Processing proteins"):
        prot_id_dict[prot_id] = {
            'cctyper': {'operon': None, 'type': None, 'subtype': None}, 
            'deffinder': {'operon': None, 'type': None, 'subtype': None}, 
            'padloc': {'operon': None, 'type': None, 'subtype': None}
        }
    
    # Fill dictionary using vectorized operations where possible
    print("[INFO] Filling protein annotations (this may take a while)...", file=sys.stderr)
    for idx, row in tqdm(merged_gff_df.iterrows(), total=len(merged_gff_df), desc="Filling annotations"):
        prot_id = row['prot_id']
        prot_id_dict[prot_id]['cctyper']['operon'] = row['operon cctyper']
        prot_id_dict[prot_id]['cctyper']['type'] = row['type cctyper']
        prot_id_dict[prot_id]['cctyper']['subtype'] = row['subtype cctyper']
        
        prot_id_dict[prot_id]['deffinder']['operon'] = row['operon deffinder']
        prot_id_dict[prot_id]['deffinder']['type'] = row['type deffinder']
        prot_id_dict[prot_id]['deffinder']['subtype'] = row['subtype deffinder']
        
        prot_id_dict[prot_id]['padloc']['operon'] = row['operon padloc']
        prot_id_dict[prot_id]['padloc']['type'] = row['type padloc']
        prot_id_dict[prot_id]['padloc']['subtype'] = row['subtype padloc']
    
    prot_id_df = pd.DataFrame.from_dict(prot_id_dict, orient='index')
    
    # Group proteins by operon
    print("[INFO] Grouping proteins by operon...", file=sys.stderr)
    operon_groups = defaultdict(set)
    for prot_id, row in tqdm(prot_id_df.iterrows(), total=len(prot_id_df), desc="Grouping operons"):
        for source in ['cctyper', 'deffinder', 'padloc']:
            operon_id = row[source]['operon']
            if pd.notnull(operon_id):
                operon_groups[operon_id].add(prot_id)
    
    operon_groups_list_dict = {key: list(value) for key, value in operon_groups.items()}
    operon_groups_df = pd.DataFrame(
        list(operon_groups_list_dict.items()),
        columns=['operon_id', 'prot_ids']
    )
    
    # Determine merged info
    print(f"[INFO] Determining merged info for {len(operon_groups_df)} operons...", file=sys.stderr)
    merged_info_df = operon_groups_df['prot_ids'].progress_apply(
        lambda x: determine_merged_info(x, prot_id_df)
    )
    
    merged_info_dict = {}
    for idx, row in operon_groups_df.iterrows():
        merged_info = merged_info_df.loc[idx]
        for prot_id in row['prot_ids']:
            merged_info_dict[prot_id] = {
                'merged_type': merged_info['merged_type'],
                'merged_subtype': merged_info['merged_subtype'],
                'merged_operon': merged_info['merged_operon']
            }
    
    # Create a new data frame with the merged information
    merged_info_df_final = pd.DataFrame.from_dict(merged_info_dict, orient='index')
    
    # Merge back
    print("[INFO] Merging unified annotations back to GFF...", file=sys.stderr)
    data_merged = merged_gff_df.merge(merged_info_df_final, left_on='prot_id', right_index=True, how="left")
    all_gff_merged = all_gff_df.merge(
        data_merged[["prot_id", "merged_type", "merged_subtype", "merged_operon"]],
        on="prot_id", how="left"
    )
    
    # Save the merged data frame
    print("[INFO] Saving unified merged file (all.merged.unified)...", file=sys.stderr)
    all_gff_merged.to_csv('all.merged.unified', index=False)
    
    print("\n" + "="*80, file=sys.stderr)
    print("[PROCESSING STAGE 4] Categorizing and summarizing systems", file=sys.stderr)
    print("="*80, file=sys.stderr)
    
    # Categorize systems
    print("[INFO] Categorizing defense systems...", file=sys.stderr)
    antidefense_types = defense_finder_df[
        defense_finder_df["activity"] == "Antidefense"
    ]["type"].unique()
    pdc_types = [n for n in padloc_df["system"].unique() if n.startswith("PDC-")]
    
    all_gff_nopdc_noanti = all_gff_merged[
        ~all_gff_merged["merged_type"].isin(list(pdc_types) + list(antidefense_types))
    ].copy()
    all_gff_df_pdc = all_gff_merged[all_gff_merged["merged_type"].isin(pdc_types)].copy()
    all_gff_df_anti = all_gff_merged[all_gff_merged["merged_type"].isin(antidefense_types)].copy()
    
    print(f"[INFO] Found {len(all_gff_nopdc_noanti[all_gff_nopdc_noanti['merged_operon'].notna()].groupby('merged_operon'))} defense systems", file=sys.stderr)
    print(f"[INFO] Found {len(all_gff_df_pdc[all_gff_df_pdc['merged_operon'].notna()].groupby('merged_operon'))} PDC systems", file=sys.stderr)
    print(f"[INFO] Found {len(all_gff_df_anti[all_gff_df_anti['merged_operon'].notna()].groupby('merged_operon'))} antidefense systems", file=sys.stderr)
    
    # Add category tags to operons
    print("[INFO] Adding category tags...", file=sys.stderr)
    all_gff_df_pdc["merged_operon"] = all_gff_df_pdc.apply(
        lambda row: row["merged_operon"] + "pdc" + str(row["seqid"])
        if pd.notnull(row["merged_operon"]) else row["merged_operon"],
        axis=1
    )
    
    all_gff_df_anti["merged_operon"] = all_gff_df_anti.apply(
        lambda row: row["merged_operon"] + "antidefense" + str(row["seqid"])
        if pd.notnull(row["merged_operon"]) else row["merged_operon"],
        axis=1
    )
    
    all_gff_nopdc_noanti["merged_operon"] = all_gff_nopdc_noanti.apply(
        lambda row: row["merged_operon"] + "defense" + str(row["seqid"])
        if pd.notnull(row["merged_operon"]) else row["merged_operon"],
        axis=1
    )
    
    # Rename columns for categorized data
    all_gff_df_pdc = all_gff_df_pdc.rename(
        columns={'merged_type': 'PDC_type', 'merged_subtype': 'PDC_subtype'}
    )
    all_gff_df_anti = all_gff_df_anti.rename(
        columns={'merged_type': 'Antidefense_type', 'merged_subtype': 'Antidefense_subtype'}
    )
    
    all_gff_merged = pd.concat([all_gff_nopdc_noanti, all_gff_df_pdc, all_gff_df_anti])
    all_gff_merged["length"] = all_gff_merged["end"] - all_gff_merged["start"]
    
    # Create summary
    print("[INFO] Creating defense systems summary...", file=sys.stderr)
    merged_data = []
    for operon, group in tqdm(all_gff_merged.groupby('merged_operon'), desc="Summarizing systems"):
        group = group.sort_values(by='start')
        start = group['start'].iloc[0]
        end = group['end'].iloc[-1]
        merged_data.append({
            'Acc': group['seqid'].iloc[0],
            'Type': group['merged_type'].iloc[0],
            'Unified_subtype_name': group['merged_subtype'].iloc[0],
            'PDC_type': group['PDC_type'].iloc[0] if 'PDC_type' in group.columns else '',
            'PDC_unified_subtype_name': group['PDC_subtype'].iloc[0] if 'PDC_subtype' in group.columns else '',
            'Antidefense_type': group['Antidefense_type'].iloc[0] if 'Antidefense_type' in group.columns else '',
            'Antidefense_unified_subtype_name': group['Antidefense_subtype'].iloc[0] if 'Antidefense_subtype' in group.columns else '',
            'Start': start,
            'End': end,
            'sys_length': end - start,
            'gene_length': sum(group['length']),
            'defense_num_genes': len(group),
        })
    
    merged_df = pd.DataFrame(merged_data)
    merged_df = merged_df.fillna("")
    
    # Categorize final results
    merged_df_nopdc_noanti = merged_df[merged_df["Type"] != ""]
    merged_df_pdc = merged_df[merged_df["PDC_type"] != ""]
    merged_df_anti = merged_df[merged_df["Antidefense_type"] != ""]
    
    # Rename columns
    merged_df_pdc = merged_df_pdc.rename(columns={
        'sys_length': 'PDC_sys_length',
        'gene_length': 'PDC_gene_length',
        'defense_num_genes': 'PDC_num_genes',
        'Start': 'PDC_start',
        'End': 'PDC_end'
    })
    
    merged_df_anti = merged_df_anti.rename(columns={
        'Antidefense_type': 'Anti_type',
        'Antidefense_unified_subtype_name': 'Anti_unified_subtype_name',
        'sys_length': 'Anti_sys_length',
        'gene_length': 'Anti_gene_length',
        'defense_num_genes': 'Anti_num_genes',
        'Start': 'Anti_start',
        'End': 'Anti_end'
    })
    
    merged_df = pd.concat([merged_df_nopdc_noanti, merged_df_pdc, merged_df_anti])
    
    print("[INFO] Saving final results (defense_results.tsv)...", file=sys.stderr)
    merged_df.to_csv('defense_results.tsv', sep="\t", index=False)
    
    print("\n" + "="*80, file=sys.stderr)
    print("[SUCCESS] Merge complete!", file=sys.stderr)
    print("="*80, file=sys.stderr)
    print(f"Total systems in final output: {len(merged_df)}", file=sys.stderr)
    print(f"  - Defense systems: {len(merged_df_nopdc_noanti)}", file=sys.stderr)
    print(f"  - PDC systems: {len(merged_df_pdc)}", file=sys.stderr)
    print(f"  - Antidefense systems: {len(merged_df_anti)}", file=sys.stderr)
    print("\nOutput files generated:", file=sys.stderr)
    print("  - all.merged", file=sys.stderr)
    print("  - all.merged.unified", file=sys.stderr)
    print("  - defense_results.tsv", file=sys.stderr)
    print("="*80, file=sys.stderr)

if __name__ == "__main__":
    main()