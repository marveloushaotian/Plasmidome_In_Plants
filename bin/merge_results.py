import pandas as pd

# Function to merge files based on Contig_ID
def merge_file(input_file, output_file):
    data = pd.read_csv(input_file, sep='\t')
    data = data.fillna('')
    # Adjust merge to skip empty values
    merged_data = data.groupby('Contig_ID').agg(lambda x: ','.join(filter(None, map(str, x)))).reset_index()
    merged_data.to_csv(output_file, sep='\t', index=False)
    print(f'Data merged and saved to {output_file}')

# Paths for defense results
input_defense = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/Defense_Results.tsv'
output_defense = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/Defense_Results_Final.tsv'

# Paths for AMR results
input_amr = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/AMR_Results.tsv'
output_amr = '~/Warehouse/GitHub/Plasmidome_In_Plants/Collect/4395/AMR_Results_Final.tsv'

# Merge both files
merge_file(input_defense, output_defense)
merge_file(input_amr, output_amr)
