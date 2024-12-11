import os
import csv
import numpy as np
import pandas as pd

# Set the main folder
main_folder_path = '/blue/salemi/share/varcovid/GenSeq'
#main_folder_path = '/Users/utente/Desktop/Varcovid/GenSeq'
# Function To load the data
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
#Datset, Dataset_w_list = load_data(main_folder_path+'/Dataset', [1])
# Calculate amino acid frequencies
print('Calculate amino acid frequencies')
all_sequences = ''.join(Datset['Sequence'])
amino_acid_counts = pd.Series(list(all_sequences)).value_counts()
total_amino_acids = amino_acid_counts.sum()
amino_acid_freqs = amino_acid_counts / total_amino_acids

# Function to generate weighted random amino acid sequences
def generate_weighted_random_sequence(start_seq):
    remaining_length = 701 - len(start_seq)
    additional_seq = np.random.choice(amino_acid_freqs.index, size=remaining_length, p=amino_acid_freqs.values)
    return start_seq + ''.join(additional_seq)

# Generate random sequences and save to CSV
unique_start = np.unique([seq[:14] for seq in Datset['Sequence']])
for start in unique_start:
    print('Generation')
    generated_sequences = [generate_weighted_random_sequence(start) for _ in range(100)]
    csv_filename = os.path.join(main_folder_path, 'Generation_Sequences', 'random_weighted', f'generated_sequences_{start}.csv')
    with open(csv_filename, mode='w', newline='', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Generated Sequence'])  # Write a header row
        for sequence in generated_sequences:
            csv_writer.writerows([[sequence]])  # Write all sequences