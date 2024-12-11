import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Path to the main folder where all models are stored
main_folder = "/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated"

# List to store dataframes from each model
df_list = []

# Columns to clean and convert to numeric
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

                # Drop unnecessary columns
                df = df.drop(columns=['gen_seq', 'E_score', 'Bit_Score'], errors='ignore')

                # Add a new column to identify which model this data belongs to
                df['Model'] = model_folder

                # Clean object columns by removing brackets and converting to numeric
                df[columns_to_convert] = df[columns_to_convert].apply(
                    lambda x: x.str.replace(r'\[|\]', '', regex=True).astype(float)
                )

                # Append the dataframe to the list
                df_list.append(df)

# Combine all dataframes into one
combined_df = pd.concat(df_list, ignore_index=True)
# Filter columns that start with 'lev_'
lev_columns = [col for col in combined_df.columns if col.startswith('lev_')]

# Select only 'lev_' columns and 'PAM30'
correlation_df = combined_df[lev_columns + ['PAM30']]

# Calculate the correlation matrix
correlation_matrix = correlation_df.corr()

# Extract correlations between 'lev_' columns and 'PAM30'
correlation_with_pam30 = correlation_matrix['PAM30'][lev_columns]

# Plot the correlations as a heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(correlation_with_pam30.to_frame(), annot=True, cmap='coolwarm', cbar=True)
plt.title("Correlation between 'lev_' columns and PAM30")
plt.xlabel("PAM30")
plt.ylabel("lev_ Columns")
plt.show()
