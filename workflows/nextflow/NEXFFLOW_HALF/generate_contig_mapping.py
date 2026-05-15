import os
import glob
import csv

input_dir = "PREPARED_FROM_TK_ALL_LAEBLED"
output_file = "contig_mapping.csv"

def generate_mapping():
    print(f"Scanning {input_dir} for FASTA files...")
    fasta_files = glob.glob(os.path.join(input_dir, "*.fasta"))
    print(f"Found {len(fasta_files)} files.")

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File_Prefix", "Contig_Name"])

        for fasta_path in fasta_files:
            filename = os.path.basename(fasta_path)
            prefix = os.path.splitext(filename)[0]
            
            try:
                with open(fasta_path, 'r') as f:
                    for line in f:
                        if line.startswith(">"):
                            contig_name = line.strip()[1:] # Remove > and newline
                            # If there are spaces, we usually take the whole line as the ID in simple mappings,
                            # or just the first word. Given the example header had no spaces, taking the whole line is safe.
                            writer.writerow([prefix, contig_name])
            except Exception as e:
                print(f"Error reading {fasta_path}: {e}")

    print(f"Mapping generated in {os.path.abspath(output_file)}")

if __name__ == "__main__":
    generate_mapping()
