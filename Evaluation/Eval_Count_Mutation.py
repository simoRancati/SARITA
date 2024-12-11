import argparse
from utils import *

# Set up command-line arguments
parser = argparse.ArgumentParser(description="Process sequences and calculate scores.")
parser.add_argument('--dir_week', type=str, required=True, help="Directory containing dataset files")
parser.add_argument('--base_path', type=str, required=True, help="Base path for generated sequences")
parser.add_argument('--output_file_path', type=str, required=True, help="Output path for saving the summary results")
parser.add_argument('--count_output_file', type=str, required=True, help="Output path for saving the amino acid count file")

args = parser.parse_args()

# Use command-line arguments
dir_week = args.dir_week
base_path = args.base_path
output_file_path = args.output_file_path
count_output_file = args.count_output_file

# Valid and target amino acids
valid_chars = set("ACDEFGHIKLMNPQRSTVWYX")
target_aminoacids = set("UO")  # Replace with your amino acids of interest

# Load the training dataset to get the first 14 amino acids
print('Downloading the Training Set')
Dataset, Dataset_w_list = load_data(dir_week, [1])  # Training set
train_step = Dataset.iloc[:, 1:len(Dataset.columns)].to_numpy()
seq_train_unique = np.unique(Dataset.iloc[:, 1:len(Dataset.columns)].to_numpy())

# Define the original Wuhan sequence
seq_origin = origin_spike()

# Extract the first unique sequences
start_seq = []
for i in range(len(train_step)):
    example = train_step[i]
    seq = example[0]
    start_seq.append(seq[:14])

unique_start = np.unique(start_seq)

# Open the count file in write mode
with open(count_output_file, 'w') as count_file:
    # Process each unique starting sequence
    for i in range(len(unique_start)):
        summary_dict = {
            "gen_seq": [],
            "len_seq": [],
            "PAM30": [],
            "BLOSUM62": [],
            "Mutation": [],
            "Counting": []
        }
        target_aminoacids_count = 0  # Initialize count of sequences containing target amino acids

        # Load the generated sequence data
        gen_seq_df = pd.read_csv(base_path + "/generated_sequences_" + str([unique_start[i]]) + ".csv")
        gen_seq = list(gen_seq_df['Generated Sequence'])

        for f in range(len(gen_seq)):
            seq_gen_filt = gen_seq[f]

            if '<unk>' in seq_gen_filt:
                seq_gen_filt = seq_gen_filt.replace('<unk>', '')  # Remove any <unk> tokens

            # Increment count if the sequence contains one of the target amino acids
            if any(amino_acid in seq_gen_filt for amino_acid in target_aminoacids):
                target_aminoacids_count += 1

            if not is_valid_sequence(seq_gen_filt, valid_chars):
                continue  # Skip invalid sequence

            summary_dict["gen_seq"].append(seq_gen_filt)
            summary_dict["len_seq"].append(len(seq_gen_filt))

            # Limit the sequence length for alignment
            if len(seq_gen_filt) > 750:
                stop = 700
            else:
                stop = len(seq_gen_filt)

            # Calculate PAM30 score and mutations
            Pam30_score, mutations_pam = calculate_pam30_scores_mutation(seq_origin[0:stop], [seq_gen_filt], 0)
            blosum62_score = calculate_blosum_scores(seq_origin[:stop], [seq_gen_filt[:stop]])

            summary_dict["PAM30"].append(Pam30_score)
            summary_dict["BLOSUM62"].append(Pam30_score)
            summary_dict["Mutation"].append(mutations_pam)
            summary_dict["Counting"].append(len(mutations_pam[0]))

        # Create summary DataFrame and save results
        summary_df = pd.DataFrame.from_dict(summary_dict)
        summary_df.to_csv(f"{output_file_path}/summary_{str(unique_start[i])}.csv", index=False)

        # Write the count of sequences containing target amino acids to the text file
        count_file.write(f"Number of sequences containing target amino acids for {unique_start[i]}: {target_aminoacids_count}\n")

        # Print the number of sequences containing the target amino acids for each set of unique sequences
        print(f"Number of sequences containing target amino acids for {unique_start[i]}: {target_aminoacids_count}")

