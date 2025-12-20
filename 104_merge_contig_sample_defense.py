import pandas as pd

# Load the data
defense_file = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/Defense_Results_Final.tsv'
amr_file = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/AMR_Results_Final.tsv'
mapping_file = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/Contig_Sample_Mapping.csv'
output_file = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/Contig_Sample_Mapping_Final.csv'

defense_data = pd.read_csv(defense_file, sep='\t')
amr_data = pd.read_csv(amr_file, sep='\t')
mapping_data = pd.read_csv(mapping_file, sep=',')

# Merge Defense with Mapping
merged_defense = pd.merge(mapping_data, defense_data, on='Contig_ID', how='left')

# Merge AMR with Mapping
merged_amr = pd.merge(mapping_data, amr_data, on='Contig_ID', how='left')

# Concatenate both merged results
final_merged = pd.concat([merged_defense, merged_amr.drop(columns=['Sample_ID'])], axis=1)

# Save the merged results
final_merged.to_csv(output_file, sep=',', index=False)

print('Final merged data saved to', output_file)

srr_file = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/SRR_sequence_info.csv'

srr_data = pd.read_csv(srr_file, sep=',')

# Merge SRR data with Mapping
final_merge = pd.merge(final_merged, srr_data, on='Sample_ID', how='left')

# Save the final merged results
final_merge.to_csv(output_file, sep=',', index=False)

print('Final merged data with SRR info saved to', output_file)
