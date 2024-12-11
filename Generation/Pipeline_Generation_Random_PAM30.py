import os
import csv
import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices
from Bio import pairwise2

# Set the main folder
main_folder_path = '/blue/salemi/share/varcovid/GenSeq'

# Function to load the data
def load_data(dir_dataset, week_range):
    week_range = [str(x) for x in week_range]
    weeks_folder = [x for x in os.listdir(dir_dataset) if x in week_range]
    df_list = []
    w_list = []

    for week in weeks_folder:
        week_path = os.path.join(dir_dataset, week)
        for root, dirs, files in os.walk(week_path):
            for file in files:
                if file.endswith(".fasta"):
                    fasta_path = os.path.join(root, file)
                    with open(fasta_path, 'r') as f:
                        header, sequence = '', ''
                        sequences = []
                        for line in f:
                            if line.startswith('>'):
                                if header:
                                    sequences.append([header, sequence])
                                header = line.strip()
                                sequence = ''
                            else:
                                sequence += line.strip()
                        if header:
                            sequences.append([header, sequence])

                        df = pd.DataFrame(sequences, columns=['Header', 'Sequence'])
                        df['Header'] = df['Header'].str.replace('>', '', regex=False)
                        df_list.append(df)
                        w_list += [week] * df.shape[0]

    return pd.concat(df_list), w_list

# Load the Training set
print('Loading Training Set')
Datset, Dataset_w_list = load_data(main_folder_path+'/dataset_nov_2023_fasta_World', [1])

# Load BLOSUM62 matrix
blosum62 = substitution_matrices.load("PAM30")

# Function to generate weighted random amino acid sequences using BLOSUM62
def generate_sequence_using_blosum62(start_seq, length=701):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    sequence = start_seq
    while len(sequence) < length:
        last_aa = sequence[-1]
        # Get a distribution for the next amino acid based on the last amino acid in the sequence
        weights = []
        for aa in amino_acids:
            pair = (last_aa, aa)
            score = blosum62[pair] if pair in blosum62 else blosum62[(aa, last_aa)]
            weights.append(np.exp(score))  # Use exponential to convert scores to probabilities
        next_aa = np.random.choice(list(amino_acids), p=np.array(weights) / np.sum(weights))
        sequence += next_aa
    return sequence

# Generate random sequences and save to CSV
unique_start = np.unique([seq[:14] for seq in Datset['Sequence']])
for start in unique_start:
    print('Generation')
    generated_sequences = [generate_sequence_using_blosum62(start) for _ in range(100)]
    print('Writing')
    csv_filename = os.path.join(main_folder_path, 'Generation_Sequences', 'random_blosum', f'generated_sequences_{start}.csv')
    with open(csv_filename, mode='w', newline='', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Generated Sequence'])  # Write a header row
        for sequence in generated_sequences:
            csv_writer.writerows([[sequence]])  # Write all sequences