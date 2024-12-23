import os
from Bio import SeqIO
import pandas as pd

def main():
    # Define input file paths
    protein_file = r"H:\DongLab\1De_novo_antibody_sequencing\Methods_and_Results\Open_pFind\anti_HA_antibody\deglycol_anti_HA_Ab_open_pFind_validation_with_template.fasta"
    pfind_file = r"H:\DongLab\1De_novo_antibody_sequencing\Methods_and_Results\Open_pFind\anti_HA_antibody\pFind-Filtered.spectra"

    # Read FASTA file and initialize dictionaries
    record_dict = {}
    record_length_dict = {}
    file_name_dict = {}

    # Parse FASTA file
    for seq in SeqIO.parse(protein_file, "fasta"):
        record_dict[seq.id] = seq.seq
        record_length_dict[seq.id] = {id: 0 for id in range(1, len(seq.seq) + 1)}
        file_name_dict[seq.id] = []

    # Read pFind data
    pfind_data = pd.read_csv(pfind_file, sep='\t', header=0)

    # Extract required columns, including File_Name
    target_columns = ['File_Name', 'Sequence', 'Modification', 'Proteins', 'Positions']
    pfind_target_data = pfind_data[target_columns]

    # Process pFind data
    for index, row in pfind_target_data.iterrows():
        File_Name, Sequence, Modification, Proteins, Positions = row
        
        if Proteins.startswith("CON_"):
            continue

        seq_length = len(Sequence)
        proteins = Proteins.split("/")
        positions_list = Positions.split("/")  # Split by slash

        for idx, one_protein in enumerate(proteins):
            if idx >= len(positions_list):  # Ensure index is not out of bounds
                break

            one_protein_sequence = record_dict.get(one_protein, None)
            if not one_protein_sequence:
                continue

            # Process corresponding Positions
            pos = positions_list[idx]
            if not pos.strip():  # Skip empty strings
                continue
            
            # Extract position number
            position_parts = pos.split(',')
            try:
                idx_sequence = int(position_parts[0])  # Try to convert the first part to an integer
            except ValueError:
                print(f"Cannot convert position '{position_parts[0]}' to integer, skipping.")
                continue

            # Update coverage count
            for sequence_position in range(idx_sequence + 1, seq_length + idx_sequence + 1):
                record_length_dict[one_protein][sequence_position] += 1

    # Create result DataFrame
    columns = ['protein', 'position', 'AA', 'spectra_number']
    df = pd.DataFrame(columns=columns)

    for protein_name, coverage_dict in record_length_dict.items():
        AA = record_dict.get(protein_name)
        for position, coverage_number in coverage_dict.items():
            if position <= len(AA):  # Ensure position is within the sequence range
                df.loc[len(df)] = [protein_name, position, AA[int(position) - 1], coverage_number]

    # Modify output file path
    output_dir = r"H:\DongLab\1De_novo_antibody_sequencing\Methods_and_Results\Open_pFind\anti_HA_antibody"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'result.xls')

    # Save results to file
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()