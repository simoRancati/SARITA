import argparse
from utils import *

# Set up command-line arguments
parser = argparse.ArgumentParser(description="Process sequences and calculate scores.")
parser.add_argument('--dir_week', type=str, required=True, help="Directory containing dataset files")
parser.add_argument('--dir_fasta', type=str, required=True, help="Directory containing test files")
parser.add_argument('--base_path', type=str, required=True, help="Base path for generated sequences")
parser.add_argument('--output_file_path', type=str, required=True, help="Output path for saving the summary results")
parser.add_argument('--count_output_file', type=str, required=True,
                    help="Output path for saving the amino acid count file")

args = parser.parse_args()

# Use command-line arguments
dir_week = args.dir_week
base_path = args.base_path
dir_fasta = args.dir_fasta
output_file_path = args.output_file_path
count_output_file = args.count_output_file

# Valid and target amino acids
valid_chars = set("ACDEFGHIKLMNPQRSTVWY")
target_aminoacids = set("UO")  # Replace with your amino acids of interest

# Load the training dataset to get the first 14 amino acids
print('Downloading the Training Set')
Dataset, Dataset_w_list = load_data(dir_week, [1])  # Training set
train_step = Dataset.iloc[:, 1:len(Dataset.columns)].to_numpy()
seq_train_unique = np.unique(Dataset.iloc[:, 1:len(Dataset.columns)].to_numpy())

print('Download the Test')
all_seq = read_protein_sequences_header(dir_fasta)

seq_train_unique = np.unique(Dataset.iloc[:, 1:len(Dataset.columns)].to_numpy())
seq_test_unique = np.unique(list(all_seq.values()))

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
        print('start the analysis for ' + str(unique_start[i]))
        summary_dict = {
            "gen_seq": [],
            "len_seq": [],
            "PAM30": [],
            "PAM30_WHU": [],
            "Mutation_Best_Gen": [],
            "Counting_Mutation_Best_Gen": [],
            "Mutation_Whuan_Best": [],
            "Counting_Mutation_Whuan_Best": [],
            "Mutation_Whuan_Gen": [],
            "Counting_Mutation_Whuan_Gen": [],
            "Mutation_intersection_best_gen": [],
            "Counting_Mutation_intersection_best_gen": [],
            'Min_dist_lev': [],
            'Seq_with_min_dist': []
        }
        target_aminoacids_count = 0  # Initialize count of sequences containing target amino acids

        # Load the generated sequence data
        gen_seq_df = pd.read_csv(base_path + "/generated_sequences_" + str([unique_start[i]]) + ".csv")
        gen_seq = list(gen_seq_df['Generated Sequence'])

        for f in range(len(gen_seq)):
            print(f)
            seq_gen_filt = gen_seq[f]

            if '<unk>' in seq_gen_filt:
                seq_gen_filt = seq_gen_filt.replace('<unk>', '')  # Remove any <unk> tokens

            # Increment count if the sequence contains one of the target amino acids
            if any(amino_acid in seq_gen_filt for amino_acid in target_aminoacids):
                target_aminoacids_count += 1

            if not is_valid_sequence(seq_gen_filt, valid_chars):
                continue  # Skip invalid sequence

            if len(seq_gen_filt) > 750:
                stop = 700
            else:
                stop = len(seq_gen_filt)

            # Compute Levenshtein distances and find the closest test sequence
            print('Compute_lev_dist')
            min_distance = float('inf')
            closest_seq_test = None
            for seq_test_real in seq_test_unique:
                distance = calculate_distance(seq_test_real[0:stop], seq_gen_filt)
                if distance < min_distance:
                    min_distance = distance
                    closest_seq_test = seq_test_real

            # Use the closest sequence for PAM30 and BLOSUM62 calculations
            summary_dict["gen_seq"].append(seq_gen_filt)
            summary_dict["len_seq"].append(len(seq_gen_filt))
            summary_dict['Min_dist_lev'].append(min_distance)
            summary_dict['Seq_with_min_dist'].append(closest_seq_test)

            # Calculate PAM30 score and mutations
            Pam30_score  = calculate_pam30_scores(closest_seq_test[:stop], [seq_gen_filt[:stop]])
            Pam30_whuan_score = calculate_pam30_scores(seq_origin[:stop], [seq_gen_filt[:stop]])

            summary_dict["PAM30"].append(Pam30_score)
            summary_dict["PAM30_WHU"].append(Pam30_whuan_score)

            # Mutation

            seq1 = closest_seq_test[0:len(seq_gen_filt)]  # The best
            seq2 = seq_gen_filt  # the generated
            seq3 = seq_origin[0:len(seq_gen_filt)]  # origin

            alignments_best_gen = pairwise2.align.globalxx(seq1, seq2)
            alignments_whuan_whuan_best = pairwise2.align.globalxx(seq3, seq1)
            alignments_whuan_whuan_gen = pairwise2.align.globalxx(seq3, seq2)

            mutations_best_gen = find_mutations(alignments_best_gen[0])
            mutations_whuan_whuan_best = find_mutations(alignments_whuan_whuan_best[0])
            mutations_whuan_whuan_gen = find_mutations(alignments_whuan_whuan_gen[0])

            summary_dict["Mutation_Best_Gen"].append(mutations_best_gen)
            summary_dict["Counting_Mutation_Best_Gen"].append(len(mutations_best_gen))
            summary_dict["Mutation_Whuan_Best"].append(mutations_whuan_whuan_best)
            summary_dict["Counting_Mutation_Whuan_Best"].append(len(mutations_whuan_whuan_best))
            summary_dict["Mutation_Whuan_Gen"].append(mutations_whuan_whuan_gen)
            summary_dict["Counting_Mutation_Whuan_Gen"].append(len(mutations_whuan_whuan_gen))

            tolerance = 2

            # Confronto delle mutazioni con tolleranza
            matched_mutations = find_matching_mutations_with_tolerance(mutations_whuan_whuan_gen,
                                                                       mutations_whuan_whuan_best, tolerance)

            summary_dict["Mutation_intersection_best_gen"].append(matched_mutations)
            summary_dict["Counting_Mutation_intersection_best_gen"].append(len(matched_mutations))

        # Create summary DataFrame and save results
        summary_df = pd.DataFrame.from_dict(summary_dict)
        summary_df.to_csv(f"{output_file_path}/summary_{str(unique_start[i])}.csv", index=False)

        # Write the count of sequences containing target amino acids to the text file
        count_file.write(
            f"Number of sequences containing target amino acids for {unique_start[i]}: {target_aminoacids_count}\n")

        # Print the number of sequences containing the target amino acids for each set of unique sequences
        print(f"Number of sequences containing target amino acids for {unique_start[i]}: {target_aminoacids_count}")