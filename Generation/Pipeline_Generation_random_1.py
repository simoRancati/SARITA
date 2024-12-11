import os
import csv
import numpy as np
import pandas as pd

# Set the main folder
main_folder_path = '/blue/salemi/share/varcovid/GenSeq'

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
Datset, Dataset_w_list = load_data(main_folder_path+'/dataset_nov_2023_fasta_World', [1])
train_step1 =Datset.iloc[:, 1:len(Datset.columns)].to_numpy()

# Find unique start
start_seq = []
for i in range(len(train_step1)):
  example = train_step1[i]
  seq = example[0]
  start_seq.append(seq[:14])

unique_start = np.unique(start_seq)
print(unique_start)

# Function to generate random amino acid sequences
def generate_random_sequence(start_seq):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # Amino acids (excluding U, B, J, Z, X)
    return start_seq + ''.join(np.random.choice(list(amino_acids), 701 - len(start_seq)))

# Generate random sequences and save to CSV
for start in unique_start:
    generated_sequences = [generate_random_sequence(start) for _ in range(100)]
    csv_filename = os.path.join(main_folder_path, 'Generation_Sequences', 'random', f'generated_sequences_{start}.csv')
    with open(csv_filename, mode='w', newline='', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Generated Sequence'])  # Write a header row
        for sequence in generated_sequences:
            csv_writer.writerows([[sequence]])  # Write all sequences