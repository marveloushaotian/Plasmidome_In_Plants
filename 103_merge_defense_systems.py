#!/usr/bin/env python3
"""
Merge defense system annotations from PADLOC, DefenseFinder, and CCTyper
"""

import pandas as pd
import sys
from collections import defaultdict
import ast

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
    """Process DefenseFinder data"""
    # Apply unified naming
    df['unified type'] = df['type'].replace(mapping_dicts['defense_finder_type'])
    df['unified subtype'] = df['subtype'].replace(mapping_dicts['defense_finder_subtype'])
    df['sys_id'] = df['sys_id'] + "|" + df["sys_beg"].astype(str).str.rsplit('_', n=1).str[0]

    # Split protein lists
    df['protein_in_syst'] = df['protein_in_syst'].str.split(',')
    df['name_of_profiles_in_sys'] = df['name_of_profiles_in_sys'].str.split(',')

    # Filter multi-contig systems
    df['contigs_in_sys'] = df['protein_in_syst'].apply(
        lambda x: len(set([item.rsplit('_',1)[0] for item in x]))
    )
    df = df[df['contigs_in_sys'] == 1]
    df.drop('contigs_in_sys', axis=1, inplace=True)

    # Expand to one row per protein
    expanded_df = df.apply(
        lambda x: pd.Series(list(zip(x['protein_in_syst'], x['name_of_profiles_in_sys']))),
        axis=1
    ).stack().reset_index(level=1, drop=True)

    expanded_df = expanded_df.apply(pd.Series)
    expanded_df.columns = ['protein_in_syst', 'name_of_profiles_in_sys']

    expanded_df = df.drop(['protein_in_syst', 'name_of_profiles_in_sys'], axis=1).merge(
        expanded_df, left_index=True, right_index=True
    )

    # Split fragmented systems
    expanded_df['seqid'] = expanded_df['sys_id'].str.split('|').str[1]
    expanded_df['prot_ord'] = expanded_df['protein_in_syst'].astype(str).str.rsplit('_',n=1).str[1].astype(int)
    expanded_df = expanded_df.sort_values(['sys_id','prot_ord']).reset_index(drop=True)

    prev_sysid = None
    prev_ord = None
    cnt = 0
    for idx, row in expanded_df.iterrows():
        if prev_sysid:
            if row["sys_id"] == prev_sysid:
                if row["prot_ord"] - prev_ord > 5:
                    prev_sysid = row["sys_id"]
                    cnt = cnt + 1
                    expanded_df.at[idx, 'sys_id'] = expanded_df.at[idx, 'sys_id'] + "|" + str(cnt)
                else:
                    if cnt > 0:
                        prev_sysid = row["sys_id"]
                        expanded_df.at[idx, 'sys_id'] = expanded_df.at[idx, 'sys_id'] + "|" + str(cnt)
            else:
                cnt = 0
                prev_sysid = row["sys_id"]
            prev_ord = row["prot_ord"]
        else:
            prev_sysid = row["sys_id"]
            prev_ord = row["prot_ord"]

    expanded_df.drop(['seqid','prot_ord'], axis=1, inplace=True)
    return expanded_df

def process_padloc(df, mapping_dicts):
    """Process PADLOC data"""
    df['unified type'] = df['system'].replace(mapping_dicts['padloc_type'])
    df['unified subtype'] = df['system'].replace(mapping_dicts['padloc_subtype'])
    df["system.number"] = df['seqid'] + "|" + df["system.number"].astype(str)
    return df

def process_cctyper(df, mapping_dicts):
    """Process CCTyper data"""
    df['unified type'] = df['Prediction'].replace(mapping_dicts['cctyper_type'])
    df['unified subtype'] = df['Prediction'].replace(mapping_dicts['cctyper_subtype'])

    # Convert string representations to lists
    df['Genes'] = df['Genes'].apply(ast.literal_eval)
    df['Prot_IDs'] = df['Prot_IDs'].apply(ast.literal_eval)

    # Expand to one row per protein
    expanded_df = df.apply(
        lambda x: pd.Series(list(zip(x['Genes'], x['Prot_IDs']))),
        axis=1
    ).stack().reset_index(level=1, drop=True)

    expanded_df = expanded_df.apply(pd.Series)
    expanded_df.columns = ['Genes', 'Prot_IDs']

    expanded_df = df.drop(['Genes', 'Prot_IDs'], axis=1).merge(
        expanded_df, left_index=True, right_index=True
    )

    return expanded_df

def aggregate_merge(main_df, merge_df, left_on, right_on, new_columns):
    """Merge with aggregation"""
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
    group_info = prot_id_df.loc[prot_id_group]

    unique_types = []
    unique_subtypes = []
    operon_names = {'cctyper': set(), 'deffinder': set(), 'padloc': set()}

    for idx, row in group_info.iterrows():
        for source in ['cctyper', 'deffinder', 'padloc']:
            if pd.notnull(row[source]['type']):
                for z in row[source]['type'].split(','):
                    if z not in unique_types:
                        unique_types.append(z)
            if pd.notnull(row[source]['subtype']):
                for z in row[source]['subtype'].split(','):
                    if z not in unique_subtypes:
                        unique_subtypes.append(z)
            if pd.notnull(row[source]['operon']):
                operon_names[source].add(row[source]['operon'])

    preference_order = ['cctyper', 'deffinder', 'padloc']
    for preference in preference_order:
        if operon_names[preference]:
            merged_operon = list(operon_names[preference])[0]
            break
    else:
        merged_operon = None

    # Remove DMS and DRT if other options available
    if len(unique_types) > 1 and ('DMS' in unique_types or 'DRT' in unique_types):
        unique_types = [t for t in unique_types if t not in ['DMS', 'DRT']]
    if len(unique_subtypes) > 1 and ('DMS' in unique_subtypes or 'DRT' in unique_subtypes):
        unique_subtypes = [s for s in unique_subtypes if s not in ['DMS', 'DRT']]

    merged_type = select_preferred_value(list(unique_types), preference_order)
    merged_subtype = select_preferred_value(list(unique_subtypes), preference_order)

    return pd.Series({
        'merged_type': merged_type,
        'merged_subtype': merged_subtype,
        'merged_operon': merged_operon
    })

def main():
    print("[INFO] Loading input files...", file=sys.stderr)

    # Load data
    all_gff_df = pd.read_csv('all_transformed_gff.tsv', sep='\t')
    padloc_df = pd.read_csv('final_padloc.tsv', sep='\t')
    defense_finder_df = pd.read_csv('final_defensefinder.tsv', sep='\t')
    cas_operons_df = pd.read_csv('final_cctyper.tsv', sep="\t")
    defense_systems_mapping = load_mapping_file('Defense Systems Name List.xlsx')

    print(f"[INFO] Loaded {len(all_gff_df)} GFF entries", file=sys.stderr)
    print(f"[INFO] Loaded {len(padloc_df)} PADLOC results", file=sys.stderr)
    print(f"[INFO] Loaded {len(defense_finder_df)} DefenseFinder results", file=sys.stderr)
    print(f"[INFO] Loaded {len(cas_operons_df)} CCTyper results", file=sys.stderr)

    # Create mappings
    mapping_dicts = create_mapping_dicts(defense_systems_mapping)

    # Process each tool's output
    print("[INFO] Processing DefenseFinder data...", file=sys.stderr)
    expanded_defense_finder_df = process_defense_finder(defense_finder_df, mapping_dicts)

    print("[INFO] Processing PADLOC data...", file=sys.stderr)
    padloc_df = process_padloc(padloc_df, mapping_dicts)

    print("[INFO] Processing CCTyper data...", file=sys.stderr)
    expanded_cas_operons_df = process_cctyper(cas_operons_df, mapping_dicts)

    # Merge all data
    print("[INFO] Merging annotations...", file=sys.stderr)
    merged_gff_df = all_gff_df.copy()

    merged_gff_df = aggregate_merge(
        merged_gff_df, expanded_defense_finder_df, 'prot_id', 'protein_in_syst',
        ['unified subtype', 'unified type', "sys_id"]
    )

    merged_gff_df = aggregate_merge(
        merged_gff_df, padloc_df, 'prot_id', 'target.name',
        ['unified subtype', 'unified type', "system.number"]
    )

    merged_gff_df = aggregate_merge(
        merged_gff_df, expanded_cas_operons_df, 'prot_id', 'Prot_IDs',
        ['unified subtype', 'unified type', "Operon"]
    )

    merged_gff_df = merged_gff_df.drop(columns=['protein_in_syst', 'target.name', 'Prot_IDs'], errors='ignore')
    merged_gff_df = merged_gff_df.dropna(
        subset=['unified type_x', 'unified subtype_x', 'sys_id',
                'unified type_y', 'unified subtype_y', 'system.number'],
        how='all'
    )

    # Aggregate by protein
    merged_gff_df = merged_gff_df.groupby('prot_id', as_index=False).agg({
        'seqid': lambda x: ','.join(map(str, x.dropna().unique())),
        'type': lambda x: ','.join(map(str, x.dropna().unique())),
        'start': 'min',
        'end': 'max',
        'score': lambda x: ','.join(map(str, x.dropna().unique())),
        'strand': lambda x: ','.join(map(str, x.dropna().unique())),
        'unified type_x': lambda x: ','.join(map(str, x.dropna().unique())),
        'unified subtype_x': lambda x: ','.join(map(str, x.dropna().unique())),
        'sys_id': lambda x: ','.join(map(str, x.dropna().unique())),
        'unified type_y': lambda x: ','.join(map(str, x.dropna().unique())),
        'unified subtype_y': lambda x: ','.join(map(str, x.dropna().unique())),
        'system.number': lambda x: ','.join(map(str, x.dropna().unique())),
        'unified type': lambda x: ','.join(map(str, x.dropna().unique())),
        'unified subtype': lambda x: ','.join(map(str, x.dropna().unique())),
        'Operon': lambda x: ','.join(map(str, x.dropna().unique())),
    })

    # Rename columns
    merged_gff_df.rename(columns={
        'unified type_x': 'type deffinder',
        'unified subtype_x': 'subtype deffinder',
        'sys_id': 'operon deffinder',
        'unified type_y': 'type padloc',
        'unified subtype_y': 'subtype padloc',
        'system.number': 'operon padloc',
        'unified type': 'type cctyper',
        'unified subtype': 'subtype cctyper',
        'Operon': 'operon cctyper'
    }, inplace=True)

    print("[INFO] Saving intermediate merged file...", file=sys.stderr)
    merged_gff_df.to_csv("all.merged", index=False, sep="\t")

    # Clean up padloc columns
    merged_gff_df['type padloc'] = merged_gff_df['type padloc'].str.split(',').str[0]
    merged_gff_df['subtype padloc'] = merged_gff_df['subtype padloc'].str.split(',').str[0]

    # Create protein ID dict
    print("[INFO] Unifying annotations...", file=sys.stderr)
    prot_id_dict = defaultdict(lambda: {
        'cctyper': {'operon': None, 'type': None, 'subtype': None},
        'deffinder': {'operon': None, 'type': None, 'subtype': None},
        'padloc': {'operon': None, 'type': None, 'subtype': None}
    })

    for idx, row in merged_gff_df.iterrows():
        prot_id_dict[row['prot_id']]['cctyper']['operon'] = row['operon cctyper']
        prot_id_dict[row['prot_id']]['cctyper']['type'] = row['type cctyper']
        prot_id_dict[row['prot_id']]['cctyper']['subtype'] = row['subtype cctyper']
        prot_id_dict[row['prot_id']]['deffinder']['operon'] = row['operon deffinder']
        prot_id_dict[row['prot_id']]['deffinder']['type'] = row['type deffinder']
        prot_id_dict[row['prot_id']]['deffinder']['subtype'] = row['subtype deffinder']
        prot_id_dict[row['prot_id']]['padloc']['operon'] = row['operon padloc']
        prot_id_dict[row['prot_id']]['padloc']['type'] = row['type padloc']
        prot_id_dict[row['prot_id']]['padloc']['subtype'] = row['subtype padloc']

    prot_id_df = pd.DataFrame.from_dict(prot_id_dict, orient='index')

    # Group proteins by operon
    operon_groups = defaultdict(set)
    for prot_id, row in prot_id_df.iterrows():
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
    merged_info_df = operon_groups_df['prot_ids'].apply(
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

    merged_info_df_final = pd.DataFrame.from_dict(merged_info_dict, orient='index')

    # Merge back to original data
    data_merged = merged_gff_df.merge(
        merged_info_df_final, left_on='prot_id', right_index=True, how="left"
    )
    all_gff_merged = all_gff_df.merge(
        data_merged[["prot_id", "merged_type", "merged_subtype", "merged_operon"]],
        on="prot_id", how="left"
    )

    print("[INFO] Saving unified merged file...", file=sys.stderr)
    all_gff_merged.to_csv('all.merged.unified', index=False)

    # Categorize systems
    antidefense_types = defense_finder_df[
        defense_finder_df["activity"] == "Antidefense"
    ]["type"].unique()
    pdc_types = [n for n in padloc_df["system"].unique() if n.startswith("PDC-")]

    all_gff_nopdc_noanti = all_gff_merged[
        ~all_gff_merged["merged_type"].isin(list(pdc_types) + list(antidefense_types))
    ]
    all_gff_df_pdc = all_gff_merged[all_gff_merged["merged_type"].isin(pdc_types)]
    all_gff_df_anti = all_gff_merged[all_gff_merged["merged_type"].isin(antidefense_types)]

    # Add category tags to operons
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
    for operon, group in all_gff_merged.groupby('merged_operon'):
        group = group.sort_values(by='start')
        start = group['start'].iloc[0]
        end = group['end'].iloc[-1]
        merged_data.append({
            'Acc': group['seqid'].iloc[0],
            'Type': group['merged_type'].iloc[0] if 'merged_type' in group.columns else "",
            'Unified_subtype_name': group['merged_subtype'].iloc[0] if 'merged_subtype' in group.columns else "",
            'PDC_type': group['PDC_type'].iloc[0] if 'PDC_type' in group.columns else "",
            'PDC_unified_subtype_name': group['PDC_subtype'].iloc[0] if 'PDC_subtype' in group.columns else "",
            'Antidefense_type': group['Antidefense_type'].iloc[0] if 'Antidefense_type' in group.columns else "",
            'Antidefense_unified_subtype_name': group['Antidefense_subtype'].iloc[0] if 'Antidefense_subtype' in group.columns else "",
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

    print("[INFO] Saving final results...", file=sys.stderr)
    merged_df.to_csv('defense_results.tsv', sep="\t", index=False)

    # Generate statistics
    with open('merge_stats.txt', 'w') as f:
        f.write(f"Total defense systems: {len(merged_df)}\\n")
        f.write(f"Defense systems: {len(merged_df_nopdc_noanti)}\\n")
        f.write(f"PDC systems: {len(merged_df_pdc)}\\n")
        f.write(f"Antidefense systems: {len(merged_df_anti)}\\n")

    print("[INFO] Merge complete!", file=sys.stderr)

if __name__ == "__main__":
    main()
