from utils import *

# Define the directory
dir_week ='/blue/salemi/share/varcovid/GenSeq/dataset_nov_2023_fasta_World'
dir_fasta = '/blue/salemi/share/varcovid/GenSeq/dataset_nov_2023_fasta_World/test_sequences.fasta'
base_path = '/blue/salemi/share/varcovid/GenSeq/Generation_Sequences/random_weighted'
output_file_path = '/blue/salemi/share/varcovid/GenSeq/Generation_Sequences/random_weighted/Results'  # define where save the performance
valid_chars = set("ACDEFGHIKLMNPQRSTVWY")

# Download Training set
print('Download the Train')
Dataset, Dataset_w_list = load_data(dir_week, [1])  # training set.
train_step = Dataset.iloc[:, 1:len(Dataset.columns)].to_numpy()

# Create Test set
print('Download the Test')
all_seq = read_protein_sequences_header(dir_fasta)

seq_train_unique = np.unique(Dataset.iloc[:, 1:len(Dataset.columns)].to_numpy())
seq_test_unique = np.unique(list(all_seq.values()))
# seq_test_unique = np.array([seq for seq in seq_all_unique if seq not in seq_train_unique])

# Download the EPI_ISL_402123 (Whuan Sequence)
seq_origin = origin_spike()

# Unique first Aminoacid
start_seq = []
summary_dict = {"gen_seq": [],
                "lev_mean": [],
                "lev_max": [],
                "lev_min": [],
                "lev_median": [],
                "lev_iqr": [],
                "lev_mean_te": [],
                "lev_max_te": [],
                "lev_min_te": [],
                "lev_median_te": [],
                "lev_iqr_te": [],
                "lev_whuan": [],
                "BLOSUM62": [],
                "PAM30": [],
                "BLOSUM80": [],
                "DAYHOFF": [],
                "E_score": [],
                "Bit_Score": []
                }
for i in range(len(train_step)):
    example = train_step[i]
    seq = example[0]
    start_seq.append(seq[:14])

unique_start = np.unique(start_seq)

# gen_seq = ['MFVFLVLLPLVSSQEFTNRTQLPSAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNYPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLSEFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGLSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGTIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVKGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEYVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGA']

print('Start Compute Value')
for i in range(len(unique_start)):
    f = 0
    # Load data
    gen_seq_df = pd.read_csv(base_path + "/generated_sequences_" + str(unique_start[i]) + ".csv")
    gen_seq = list(gen_seq_df['Generated Sequence'])

    distance_lev_tr = []
    distance_lev_te = []
    for f in range(len(gen_seq)):
        seq_gen_filt = gen_seq[f]
        # seq_gen_filt = gen_seq[f][:14] + gen_seq[f][19:]

        if '<unk>' in gen_seq[f]:
            seq_gen_filt = gen_seq[f].replace('<unk>', '')  # if <unk> is present

        if not is_valid_sequence(seq_gen_filt,valid_chars):
            continue  # drop this sequence

        # TRAIN
        for l in range(len(seq_train_unique)):
            seq_real = seq_train_unique[l]
            distance_lev_tr.append(calculate_distance(seq_real[0:len(seq_gen_filt)], seq_gen_filt))
        # TEST
        for g in range(len(seq_test_unique)):
            seq_test_real = seq_test_unique[g]
            distance_lev_te.append(calculate_distance(seq_test_real[0:len(seq_gen_filt)], seq_gen_filt))

        lev_dict = statistics(distance_lev_tr)
        lev_dict_test = statistics(distance_lev_te)
        dict_bio_index = Escore_bitscore(seq_origin[0:len(seq_gen_filt)], seq_gen_filt)

        # TRAIN
        summary_dict["gen_seq"].append(seq_gen_filt)
        summary_dict["lev_mean"].append(lev_dict['mean'])
        summary_dict["lev_max"].append(lev_dict['max'])
        summary_dict["lev_min"].append(lev_dict['min'])
        summary_dict["lev_iqr"].append(lev_dict['iqr'])

        # TEST
        summary_dict["lev_median"].append(lev_dict['median'])
        summary_dict["lev_mean_te"].append(lev_dict_test['mean'])
        summary_dict["lev_max_te"].append(lev_dict_test['max'])
        summary_dict["lev_min_te"].append(lev_dict_test['min'])
        summary_dict["lev_iqr_te"].append(lev_dict_test['iqr'])
        summary_dict["lev_median_te"].append(lev_dict_test['median'])

        # BIOLOGY
        summary_dict["lev_whuan"].append(calculate_distance(seq_origin[0:len(seq_gen_filt)], seq_gen_filt))
        summary_dict["BLOSUM62"].append(calculate_blosum_scores(seq_origin[0:len(seq_gen_filt)], [seq_gen_filt]))
        summary_dict["PAM30"].append(calculate_pam30_scores(seq_origin[0:len(seq_gen_filt)], [seq_gen_filt]))
        summary_dict["BLOSUM80"].append(calculate_blosum80_scores(seq_origin[0:len(seq_gen_filt)], [seq_gen_filt]))
        summary_dict["DAYHOFF"].append(calculate_dayhoff_scores(seq_origin[0:len(seq_gen_filt)], [seq_gen_filt]))
        summary_dict["E_score"].append(dict_bio_index['escore'])
        summary_dict['Bit_Score'].append(dict_bio_index['bitscore'])

    summary_df = pd.DataFrame.from_dict(summary_dict)
    summary_df.to_csv(output_file_path + "/summary_" + str(unique_start[i]) + ".csv")