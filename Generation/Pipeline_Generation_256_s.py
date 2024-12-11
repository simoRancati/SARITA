import os
import csv
import numpy as np
import pandas as pd
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer

# Set the main folder
main_folder_path = '/blue/salemi/share/varcovid/GenSeq'

# Load model and tokenizer
model = AutoModelForCausalLM.from_pretrained(main_folder_path+'/Model_FineTuned/RITA_FineTuning_runclm_256_s/checkpoint-100130', trust_remote_code=True)
tokenizer = AutoTokenizer.from_pretrained(main_folder_path+'/Model_FineTuned/RITA_FineTuning_runclm_256_s/checkpoint-100130')

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
                        # Qui leggiamo il file FASTA e lo convertiamo in un formato DataFrame
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
                        # Aggiungi l'ultima sequenza
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

# Check for GPU availability and move the model to GPU
device = "cuda" if torch.cuda.is_available() else "cpu"
model = model.to(device)

for i in range(len(unique_start)):
    # Prepare model inputs
    model_inputs = tokenizer([unique_start[i]], return_tensors="pt")
    model_inputs = {k: v.to(device) for k, v in model_inputs.items()}

    # Generate predictions using the model
    generated_ids = model.generate(**model_inputs, min_length=701, max_length=701,
                                   do_sample=True, top_k=950, repetition_penalty=1.2,
                                   num_return_sequences=100, eos_token_id=2, truncation=True)

    # Decode and print outputs
    generated_sequences = []
    for f in range(len(generated_ids)):
        sequence = tokenizer.decode(generated_ids[f], skip_special_tokens=True).replace(' ', '')
        generated_sequences.append(sequence)

    # Save all generated sequences to a CSV file
    csv_filename = main_folder_path + '/Generation_Sequences/run_clm_256_s/' + 'generated_sequences_' + str([unique_start[i]]) + '.csv'
    with open(csv_filename, mode='w', newline='', encoding='utf-8') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(['Generated Sequence'])  # Write a header row
        for f in range(len(generated_sequences)):
            csv_writer.writerows([[generated_sequences[f]]])  # Write all sequences