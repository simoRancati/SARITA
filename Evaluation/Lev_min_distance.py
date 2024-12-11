import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Path to the main folder where all models are stored
main_folder = "/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated"

# List to store dataframes from each model
df_list = []

columns_to_convert = ['BLOSUM62', 'PAM30', 'BLOSUM80', 'DAYHOFF']

# Loop through each model's folder
for model_folder in os.listdir(main_folder):
    model_path = os.path.join(main_folder, model_folder)
    if os.path.isdir(model_path):  # Check if it's a folder
        results_path = os.path.join(model_path, "Results")

        # Read all CSV files in the Results folder of the current model
        for file in os.listdir(results_path):
            if file.endswith(".csv"):
                file_path = os.path.join(results_path, file)

                # Read the CSV file
                df = pd.read_csv(file_path)

                # Drop the gen_seq column and the last two columns
                df = df.drop(columns=[df.columns[-1], df.columns[-2]])

                # Add a new column to identify which model this data belongs to
                df['Model'] = model_folder

                # Clean object columns by removing brackets and converting to numeric
                df = df.apply(lambda x: x.str.replace(r'\[|\]', '', regex=True) if x.dtype == 'object' else x)

                # Convert all object columns to numeric, after cleaning
                df[columns_to_convert] = df[columns_to_convert].apply(pd.to_numeric, errors='coerce')

                # Append the dataframe to the list
                df_list.append(df)

# Combine all dataframes into one
combined_df = pd.concat(df_list, ignore_index=True)

# Drop duplicate rows based on the 'gen_seq' and 'Model' column
combined_df = combined_df.drop_duplicates(subset=['gen_seq', 'Model'])

combined_df = combined_df.drop(columns=['gen_seq', combined_df.columns[0]])

model_mapping = {
    "SpikeGTP": "SpikeGPT2",
    "run_clm_256_s": "SARITA S",
    "rita_s": "RITA S",
    "rita_xl": "RITA XL",
    "random_weighted": "Rand Fr",
    "run_clm_256_m": "SARITA M",
    "run_clm_256_l": "SARITA L",
    "random_blosum": "Rand PAM30",
    "rita_m": "RITA M",
    "rita_l": "RITA L",
    "random": "Rand",
    "run_clm_256_xl": "SARITA XL"
}

# DataFrame
combined_df['Model'] = combined_df['Model'].map(model_mapping)

combined_df = combined_df[['lev_min_te', 'Model']]

# Funzione per calcolare la percentuale cumulativa
def cumulative_percentage(group):
    # Ordinamento dei dati per 'lev_min_te' in modo crescente
    sorted_group = group.sort_values(by='lev_min_te')
    # Calcolo della frequenza cumulativa
    sorted_group['cumulative_count'] = np.arange(1, len(sorted_group) + 1)
    # Calcolo della percentuale rispetto al totale delle sequenze nel modello
    sorted_group['percentage'] = (sorted_group['cumulative_count'] / len(group)) * 100
    return sorted_group

# Applicazione della funzione per ogni modello
cumulative_df = combined_df.groupby('Model').apply(cumulative_percentage)

# Reset dell'indice per facilitare il plotting
cumulative_df = cumulative_df.reset_index(drop=True)

# Creazione dei grafici
models = cumulative_df['Model'].unique()

# Creazione del grafico unico con legenda per ogni modello
plt.figure(figsize=(10, 6))  # Set the figure size for better visibility
for model in models:
    subset = cumulative_df[cumulative_df['Model'] == model]
    plt.plot(subset['lev_min_te'], subset['percentage'], marker='o', linestyle='-', label=model)
plt.title('Cumulative Percentage of Valid Sequences by Model')
plt.xlabel('Levenshtein Distance')
plt.ylabel('Cumulative Percentage (%)')
plt.legend(title='Model')  # Adding a legend with title
plt.grid(True)
plt.show()

import seaborn as sns

model_colors = {
    "SpikeGPT2": "#7f7f7f",
    "SARITA S": "#bcbd22",
    "RITA S": "#d62728",
    "RITA XL": "#e377c2",
    "Rand Fr": "#ff7f0e",
    "SARITA M": "#17becf",
    "SARITA L": "#ff9896",
    "Rand PAM30": "#2ca02c",
    "RITA M": "#9467bd",
    "RITA L": "#8c564b",
    "Rand": "#1f77b4",
    "SARITA XL": "#98df8a"
}

# Impostazione dello stile per un aspetto più pulito e moderno
sns.set(style="whitegrid")

# Creazione della figura con legenda per ogni modello
plt.figure(figsize=(12, 8))  # Aumenta le dimensioni per una migliore leggibilità

for model in models:
    subset = cumulative_df[cumulative_df['Model'] == model]
    plt.plot(subset['lev_min_te'], subset['percentage'], marker='o', linestyle='-', color=model_colors[model], label=model)

plt.title('Cumulative Percentage of Valid Sequences by Model', fontsize=16)
plt.xlabel('Levenshtein Distance', fontsize=14)
plt.ylabel('Cumulative Percentage (%)', fontsize=14)
plt.legend(title='Model', title_fontsize='13', fontsize='12')

# Rimozione dei bordi del grafico
ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()  # Aggiustamento automatico dei sottopannelli per riempire l'area del grafico
plt.show()
