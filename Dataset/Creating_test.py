from Bio import SeqIO

def read_fasta(file_path):
    """Read a FASTA file and return a list of sequences."""
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def write_fasta(sequences, output_file):
    """Write a list of sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq{i+1}\n{seq}\n")

def get_unique_sequences(sequences):
    """Return a set of unique sequences."""
    return set(sequences)

def get_test_sequences(all_sequences, train_sequences):
    """Return sequences that are in all_sequences but not in train_sequences."""
    train_sequences_set = set(train_sequences)
    return [seq for seq in all_sequences if seq not in train_sequences_set]

# Paths to your input FASTA files
train_fasta_path = '/blue/salemi/share/varcovid/GenSeq/dataset_nov_2023_fasta_World/1/sequences.fasta'
all_fasta_path = '/blue/salemi/share/varcovid/DeepAutoCov/spikeprot1105/spikes.fasta'
output_train_fasta_path = '/blue/salemi/share/varcovid/GenSeq/dataset_nov_2023_fasta_World/1/sequences_unique_train.fasta'
output_test_fasta_path = '/blue/salemi/share/varcovid/GenSeq/dataset_nov_2023_fasta_World/1/test_sequences.fasta'

# Read the train and all sequences
print('Read Sequences Fasta Train')
train_sequences = read_fasta(train_fasta_path)
print('Read All Sequences')
all_sequences = read_fasta(all_fasta_path)

# Get unique train sequences
unique_train_sequences = get_unique_sequences(train_sequences)

# Write the test sequences to a new FASTA file
write_fasta(unique_train_sequences, output_train_fasta_path)

print(f"Train sequences written to {output_train_fasta_path}")
unique_all_sequences = get_unique_sequences(all_sequences)
# Get test sequences (unique sequences not in train)
test_sequences = get_test_sequences(unique_all_sequences, unique_train_sequences)

# Write the test sequences to a new FASTA file
write_fasta(test_sequences, output_test_fasta_path)

print(f"Test sequences written to {output_test_fasta_path}")