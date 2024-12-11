import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import mannwhitneyu
from itertools import combinations

# Path to the main folder where all models are stored
main_folder = "/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated"

# List to store dataframes from each model
df_list = []
df_seq_list = []

#columns_to_convert = ['BLOSUM62', 'PAM30', 'BLOSUM80', 'DAYHOFF']
columns_to_convert = ['len_seq','PAM30', 'PAM30_WHU','Counting_Mutation_Best_Gen', 'Counting_Mutation_Whuan_Best','Counting_Mutation_Whuan_Gen','Counting_Mutation_intersection_best_gen']

# Loop through each model's folder
for model_folder in os.listdir(main_folder):
    model_path = os.path.join(main_folder, model_folder)
    if os.path.isdir(model_path):  # Check if it's a folder
        results_path = os.path.join(model_path, "Results_mut")

        # Read all CSV files in the Results folder of the current model
        for file in os.listdir(results_path):
            if file.endswith(".csv"):
                file_path = os.path.join(results_path, file)

                try:
                    # Read the CSV file
                    df = pd.read_csv(file_path)

                    # Drop the gen_seq column and the last column, if they exist
                    columns_to_drop = [col for col in ['gen_seq', df.columns[-1]] if col in df.columns]
                    df_seq = df[columns_to_drop]
                    df = df.drop(columns=columns_to_drop)

                    # Add a new column to identify which model this data belongs to
                    df['Model'] = model_folder
                    df_seq['Model'] = model_folder

                    # Clean object columns by removing brackets and converting to numeric
                    df = df.apply(lambda x: x.str.replace(r'\[|\]', '', regex=True) if x.dtype == 'object' else x)

                    # Convert all object columns to numeric, after cleaning
                    df[columns_to_convert] = df[columns_to_convert].apply(pd.to_numeric, errors='coerce')

                    # Append the dataframe to the list
                    df_list.append(df)
                    df_seq_list.append((df_seq))

                except KeyError as e:
                    print(f"Skipping file {file_path} due to missing columns: {e}")
                    continue

# Combine all dataframes into one
combined_df = pd.concat(df_list, ignore_index=True)
combined_df_seq = pd.concat(df_seq_list, ignore_index=True)

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

# Definisci la matrice PAM30
pam30 = substitution_matrices.load("PAM30")
def calculate_pam30_scores(reference, sequences):
    # Crea un oggetto PairwiseAligner
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = pam30

    # Lista per memorizzare i punteggi
    scores = []

    # Allineamento globale per ogni sequenza e memorizza il punteggio
    for sequence in sequences:
        # Prendi il segmento della sequenza di riferimento della stessa lunghezza della sequenza generata
        sub_reference = reference[:len(sequence)]
        score = aligner.score(sub_reference, sequence)
        scores.append(score)

    # Restituisci la lista dei punteggi
    return np.array(scores)

# Calcola il punteggio PAM30 per ogni riga del dataframe
combined_df_seq['pam30_score'] = combined_df_seq.apply(lambda row: calculate_pam30_scores(row['Seq_with_min_dist'], [row['gen_seq']])[0], axis=1)

# Calcola il punteggio mediano per ogni modello
median_scores = combined_df_seq.groupby('Model')['pam30_score'].median().reset_index()

# Mostra i risultati
print(median_scores)

combined_df['False Mutation Rate'] = (combined_df['Counting_Mutation_Whuan_Gen'] - combined_df['Counting_Mutation_intersection_best_gen']) / combined_df['Counting_Mutation_Whuan_Gen']
cols_to_convert = ['len_seq', 'PAM30', 'PAM30_WHU', 'Counting_Mutation_Best_Gen', 'Counting_Mutation_Whuan_Best', 'Counting_Mutation_Whuan_Gen', 'Counting_Mutation_intersection_best_gen', 'False Mutation Rate']
combined_df[cols_to_convert] = combined_df[cols_to_convert].apply(pd.to_numeric, errors='coerce')# Drop duplicate rows based on the 'gen_seq' and 'Model' column
#combined_df = combined_df.drop_duplicates(subset=['gen_seq', 'Model'])
# Add a new column for False Positive Rate

#combined_df = combined_df.drop(columns=['gen_seq', combined_df.columns[0]])

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

combined_df['Model'] = combined_df['Model'].map(model_mapping)
# Only select numeric columns for plotting
numeric_columns = combined_df.select_dtypes(include=['number']).columns

# Calculate the percentage of len_seq=701 for each model
len_seq_range_percentage = combined_df.groupby('Model').apply(
    lambda x: ((x['len_seq'] >= 686-30) & (x['len_seq'] <= 686+30)).mean() * 100
).reset_index(name='Len_Seq_Range_686_710_Percentage')

columns_to_convert = ['PAM30_WHU','False Mutation Rate']
# Calculate median values for each metric by model
medians = combined_df.groupby('Model')[columns_to_convert].median()
medians.to_csv("/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated/model_medians.csv")

# Define the palette once
#palette = sns.color_palette("Set2")

normalize = False  # Set to True if you want to normalize the data

# Define the new order of models
model_order = ["Rand", "Rand Fr", "Rand PAM30", "RITA S", "RITA M", "RITA L", "RITA XL",
               "SpikeGPT2", "SARITA S", "SARITA M", "SARITA L", "SARITA XL"]

# Define the specific colors
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
          '#7f7f7f', '#bcbd22', '#17becf', '#ff9896', '#98df8a']

# Perform pairwise Mann-Whitney U test for each metric and for each pair of models
for metric in numeric_columns:
    # model_pairs = list(combinations(combined_df['Model'].unique(), 2))  # Get all pairs of models
    #
    # print(f"\nMann-Whitney U test for {metric}:\n")
    #
    # for model1, model2 in model_pairs:
    #     # Get data for the two models
    #     data1 = combined_df[combined_df['Model'] == model1][metric].dropna()
    #     data2 = combined_df[combined_df['Model'] == model2][metric].dropna()
    #
    #     # Perform the Mann-Whitney U test
    #     stat, p_value = mannwhitneyu(data1, data2)
    #
    #     # Print the results
    #     print(f"Comparison: {model1} vs {model2}")
    #     print(f"U Statistic: {stat}, p-value: {p_value}\n")

    # Plot boxplots for each metric with enhanced graphics
    plt.figure(figsize=(12, 8))
    sns.boxplot(x="Model", y=metric, data=combined_df, order=model_order, palette=colors)
    plt.title(f"Comparison of {metric} Across Models", fontsize=16)
    plt.xlabel("Model", fontsize=14)
    plt.ylabel(f"{metric} (Normalized)" if normalize else metric, fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    # Removing borders
    for spine in plt.gca().spines.values():
        spine.set_visible(False)
    plt.tight_layout()
    plt.savefig(str(metric)+'.png')
    plt.show()
    print()

# NORMALIZED

# Optional Normalization (Enable one of the following)
normalize = True  # Set to True if you want to normalize the data

if normalize:
    scaler = MinMaxScaler()  # You can also use StandardScaler() for Z-score normalization
    # Select columns to normalize (all columns except 'Model')
    columns_to_normalize = combined_df.columns.drop('Model')

    # Apply the scaler only to the selected columns
    combined_df[columns_to_normalize] = scaler.fit_transform(combined_df[columns_to_normalize])

# Improved Visualization
sns.set(style="whitegrid")
palette = sns.color_palette("Set2")

# Plot boxplots for each numeric column with enhanced graphics
for metric in columns_to_convert:
    plt.figure(figsize=(12, 8))
    sns.boxplot(x="Model", y=metric, data=combined_df, palette=palette)
    plt.title(f"Comparison of {metric} Across Models", fontsize=16)
    plt.xlabel("Model", fontsize=14)
    plt.ylabel(f"{metric} (Normalized)" if normalize else metric, fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

print('Done')

