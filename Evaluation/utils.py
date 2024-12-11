import os
import pandas as pd
from Levenshtein import distance as levenshtein_distance
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio import pairwise2
import math
from collections import defaultdict

def is_valid_sequence(sequence, valid_chars):
    """Verifica se la sequenza contiene solo caratteri validi."""
    return all(char in valid_chars for char in sequence)

def Escore_bitscore(ref_seq, gen_seq):
    """
    Calculate the bit score and E-score for sequence alignments between
    a reference sequence (ref_seq) and a generated sequence (gen_seq)
    using the BLOSUM62 substitution matrix.

    This function performs a global alignment of two sequences and computes
    the bit score based on the formula from Altschul et al. (1997). The E-score
    is then calculated to estimate the number of random alignments expected to
    achieve a similar or better score by chance.

    Parameters:
    - ref_seq (str): The reference sequence to be aligned.
    - gen_seq (str): The generated sequence to be aligned.

    Returns:
    - dict: A dictionary containing the 'bitscore' and 'escore' as keys with
            their respective calculated values.

    Reference:
    Altschul, Stephen F., et al. "Gapped BLAST and PSI-BLAST: a new generation
    of protein database search programs." Nucleic acids research 25.17 (1997): 3389-3402.

    Bit Score Explanation:
    ----------------------
    The bit score (S') is a normalized score that allows for the comparison of
    alignment scores from different searches, even those using different scoring
    matrices. It is calculated using the following formula:

        bit_score = (λ * raw_score - ln(K)) / ln(2)

    where:
    - λ (lambda) is a statistical parameter dependent on the scoring system.
    - raw_score is the score of the alignment (S).
    - K is a constant used in the calculation of the bit score.

    The bit score provides a measure of the similarity between two sequences.
    Higher bit scores indicate more significant alignments. The normalization
    allows for comparison across different searches.

    Specific ranges and interpretation:
    - **High bit score**:
        * Indicates a highly significant alignment.
        * Higher values mean better alignments and greater similarity.

    - **Low bit score**:
        * Indicates a less significant alignment.
        * Lower values mean weaker alignments and less similarity.

    E-score Explanation:
    --------------------
    The E-score (or E-value) represents the number of alignments with a score
    equal to or better than the observed score that are expected to occur by
    chance alone. It provides a measure of the statistical significance of
    the alignment. The interpretation of E-scores can be summarized as follows:

    - **Small E-score (close to zero)**:
        * Example: 1e-50, 1e-100
        * Significance: Extremely significant
        * Interpretation: Highly unlikely to occur by chance. Indicates a strong
                          biological relevance.

    - **Intermediate E-score**:
        * Example: 0.001, 1
        * Significance: Moderately significant
        * Interpretation: Could be biologically relevant, but further investigation
                          is needed. An E-score of 1 means that one such alignment
                          is expected to occur by chance in the database.

    - **Large E-score (greater than 1)**:
        * Example: 10, 100
        * Significance: Not significant
        * Interpretation: Likely to occur by chance. Indicates that the alignment
                          is probably not biologically relevant.

    Specific examples:
    - **E-score of 1e-100**: Extremely significant. Indicates a strong biological relevance.
    - **E-score of 1e-20**: Very significant.
    - **E-score of 0.001**: Significant. Indicates a small chance of occurring by chance.
    - **E-score of 1**: Moderately significant. Needs further investigation.
    - **E-score of 10 or higher**: Not significant. Likely to be a random alignment.
    """

    matrix = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.globaldx(ref_seq, gen_seq, matrix)

    best_alignment = alignments[0]
    score = best_alignment[2]

    lambda_ = 0.3176
    K = 0.134

    bit_score = (lambda_ * score - math.log(K)) / math.log(2)

    m = len(ref_seq)
    n = len(gen_seq)
    E_score = m * n * math.pow(2, -bit_score)

    return {'bitscore' : bit_score,
            'escore' : E_score}

def calculate_distance(seq1, seq2):
    """
    Calculate the Levenshtein distance between two amino acid sequences.

    Args:
    seq1 (str): First amino acid sequence. label
    seq2 (str): Second amino acid sequence. output

    Returns:
    int: Levenshtein distance between seq1 and seq2.
    """
    return levenshtein_distance(seq1, seq2)

def calculate_blosum_scores(reference, sequences):
    """
    Calculate BLOSUM62 alignment scores between a reference sequence and a list of sequences.

    Input:
    - reference: A string representing the reference protein sequence.
    - sequences: A list of strings representing the protein sequences to be compared against the reference.

    Output:
    - scores: A list of floats representing the alignment scores for each sequence in the 'sequences' list,
              maintaining the order of the input sequences.
    """
    # Load the BLOSUM62 matrix
    matrix = substitution_matrices.load("BLOSUM62")

    # Create a PairwiseAligner object
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix

    # List to store the scores
    scores = []

    # Perform global alignment for each sequence and store the score
    for sequence in sequences:
        score = aligner.score(reference, sequence)
        scores.append(score)

    # Return the list of scores
    return np.array(scores)

def calculate_pam30_scores(reference, sequences):
    """
    Calculate PAM30 alignment scores between a reference sequence and a list of sequences.

    Input:
    - reference: A string representing the reference protein sequence.
    - sequences: A list of strings representing the protein sequences to be compared against the reference.

    Output:
    - scores: A list of floats representing the alignment scores for each sequence in the 'sequences' list,
              maintaining the order of the input sequences.
    """
    # Load the PAM30 matrix
    matrix = substitution_matrices.load("PAM30")

    # Create a PairwiseAligner object
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix

    # List to store the scores
    scores = []

    # Perform global alignment for each sequence and store the score
    for sequence in sequences:
        score = aligner.score(reference, sequence)
        scores.append(score)

    # Return the list of scores
    return np.array(scores)

def calculate_blosum80_scores(reference, sequences):
    """
    Calculate PAM30 alignment scores between a reference sequence and a list of sequences.

    Input:
    - reference: A string representing the reference protein sequence.
    - sequences: A list of strings representing the protein sequences to be compared against the reference.

    Output:
    - scores: A list of floats representing the alignment scores for each sequence in the 'sequences' list,
              maintaining the order of the input sequences.
    """
    # Load the BLOSUM80 matrix
    matrix = substitution_matrices.load("BLOSUM80")

    # Create a PairwiseAligner object
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix

    # List to store the scores
    scores = []

    # Perform global alignment for each sequence and store the score
    for sequence in sequences:
        score = aligner.score(reference, sequence)
        scores.append(score)

    # Return the list of scores
    return np.array(scores)

def calculate_dayhoff_scores(reference, sequences):
    """
    Calculate PAM30 alignment scores between a reference sequence and a list of sequences.

    Input:
    - reference: A string representing the reference protein sequence.
    - sequences: A list of strings representing the protein sequences to be compared against the reference.

    Output:
    - scores: A list of floats representing the alignment scores for each sequence in the 'sequences' list,
              maintaining the order of the input sequences.
    """
    # Load the DAYHOFF matrix
    matrix = substitution_matrices.load('DAYHOFF')

    # Create a PairwiseAligner object
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix

    # List to store the scores
    scores = []

    # Perform global alignment for each sequence and store the score
    for sequence in sequences:
        score = aligner.score(reference, sequence)
        scores.append(score)

    # Return the list of scores
    return np.array(scores)






def read_protein_sequences_header(file):
    """
    This function reads a FASTA file, extracts protein sequences, and their respective headers,
    and stores them in a dictionary with headers as keys and sequences as values.

    Parameters:
    - file: The path to the FASTA file.

    Returns:
    - sequences_dict: A dictionary where each key is a header and the corresponding value is the protein sequence.
    """

    sequences = []  # Initialize an empty list to store protein sequences.
    headers = []    # Initialize an empty list to store headers of the sequences.
    sequences_dict = {}  # Initialize an empty dictionary to store the sequences with headers as keys.

    with open(file, 'r', encoding='utf-8', errors='ignore') as f:
        current_sequence = ''  # A variable to hold the current sequence.

        for line in f:
            line = line.strip()  # Strip whitespace from the line.

            if line.startswith('>'):  # Check for header line.
                if current_sequence:  # If there's a current sequence, store it before starting a new one.
                    sequences.append(current_sequence)
                    headers.append(current_header)  # Store the current header
                    sequences_dict[current_header] = current_sequence  # Map the current header to its sequence
                    current_sequence = ''

                # Update the current header, excluding the first character ('>')
                current_header = line[1:]

            else:
                current_sequence += line  # Append non-header lines to the current sequence.

        if current_sequence:  # After the last line, add the final sequence.
            sequences.append(current_sequence)
            headers.append(current_header)
            sequences_dict[current_header] = current_sequence  # Map the last header to its sequence

    return sequences_dict

def load_data(dir_dataset, week_range):
    """
    Load and combine FASTA files from specified weekly data folders into a single DataFrame.

    This function processes a directory containing subdirectories for different weeks, each
    potentially holding FASTA files. It identifies the relevant week folders based on the input
    `week_range`, reads the sequences from the FASTA files within these folders, and compiles
    them into a DataFrame. Each sequence is paired with its header, and an accompanying list
    tracks the week each sequence belongs to.

    Parameters:
    - dir_dataset (str): The directory path containing subdirectories for each week.
    - week_range (list): A list of week identifiers (integers or strings) specifying which
                         weeks to include in the data processing.

    Returns:
    - pd.DataFrame: A DataFrame containing all sequences and their headers from the specified
                    weeks. Columns are 'Header' and 'Sequence'.
    - list: A list of week identifiers corresponding to each sequence in the DataFrame.

    The function filters the subdirectories in `dir_dataset` to match the entries in
    `week_range`, ensuring only the desired weeks are processed. It navigates through
    each relevant directory, identifying `.fasta` files and reading the sequence data.
    For each FASTA file, it reads and processes the headers and sequences, accumulating
    them into a DataFrame. The headers are cleaned by removing the '>' character.
    A list with week identifiers is also maintained, matching each sequence to its source week.

    Example usage:
    --------------
    dir_dataset = "path/to/dataset"
    week_range = [1, 2, 3]
    data, weeks = load_data(dir_dataset, week_range)
    """
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

def origin_spike():
    """
        Return Spike Protein isolated in Whuan.
    """
    seq_orig = 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT'


    return seq_orig

def statistics(measures):
    """
    Calculate and return common statistical measures for a given dataset.

    This function computes various statistical metrics including the maximum, minimum,
    median, mean, and interquartile range (IQR) for an array of numerical data provided
    in `measures`. These metrics help in understanding the distribution and spread of
    the data.

    Parameters:
    - measures (np.array or list): A numeric dataset for which statistics will be calculated.

    Returns:
    - dict: A dictionary containing the computed statistical values:
        'max' : Maximum value in the dataset.
        'min' : Minimum value in the dataset.
        'median' : Median value of the dataset.
        'mean' : Mean (average) value of the dataset.
        'iqr' : Interquartile range, representing the range between the 25th and 75th percentiles.

    The function leverages numpy functions to efficiently compute these statistics. It
    is suitable for analysis in scientific research, data analysis, and more contexts
    where basic descriptive statistics are needed.

    Example usage:
    --------------
    data = [10, 20, 30, 40, 50]
    stats = statistics(data)
    """
    max_value = np.max(measures)
    min_value = np.min(measures)
    median_value = np.median(measures)
    mean_value = np.mean(measures)
    q1,q3 = np.percentile(measures, [25,75])
    iqr = q3 - q1
    return {'max' : max_value,
            'min' : min_value,
            'median' : median_value,
            'mean' : mean_value,
            'iqr' : iqr}


def calculate_pam30_scores_mutation(reference, sequences, ref_start):
    """
    Calculate PAM30 alignment scores and mutations between a reference sequence and a list of sequences.

    Input:
    - reference: A string representing the reference protein sequence.
    - sequences: A list of strings representing the protein sequences to be compared against the reference.

    Output:
    - scores: A list of floats representing the alignment scores for each sequence in the 'sequences' list.
    - mutations: A list of lists, where each inner list contains tuples of (position, ref_aa, seq_aa) representing
                 the mutations found in each sequence compared to the reference.
    """

    REFSEQ_START = ref_start  # START S1

    def check_alignment_substitutions(x):
        subs = filter(lambda x: x != '',
                      map(lambda i: (''.join([x.target[i], str(REFSEQ_START + i + 1), x.query[i]])
                                     if x.target[i] != x.query[i]
                                     else ''), range(len(x.query))))
        return list(subs)
    # Load the PAM30 matrix
    matrix = substitution_matrices.load("PAM30")


    # Create a PairwiseAligner object
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix

    # Lists to store the scores and mutations
    scores = []
    mutations = []

    # Perform global alignment for each sequence
    for sequence in sequences:
        all_substitutions_predict_dict = defaultdict(int)
        # Calculate the alignment score
        score = aligner.score(reference, sequence)
        scores.append(score)

        # Retrieve the alignment details to identify mutations
        sequences_mutations = []
        alignments = aligner.align(reference, sequence)
        for k in check_alignment_substitutions(alignments[0]):
            all_substitutions_predict_dict[k] += 1

        all_substitutions = set(all_substitutions_predict_dict.keys())
        mutations.append(all_substitutions)

    # Return the list of scores and mutations
    return np.array(scores), mutations


def find_mutations(alignment):
    ref_seq, sample_seq, _, _, _ = alignment
    mutations = []

    for i, (ref_res, sample_res) in enumerate(zip(ref_seq, sample_seq), start=1):
        if ref_res != sample_res:
            if ref_res == "-":
                mutations.append(f"Inserzione {sample_res} in {i}")
            elif sample_res == "-":
                mutations.append(f"Delezione {ref_res} in {i}")
            else:
                mutations.append(f"{ref_res}{i}{sample_res}")

    return mutations

def find_matching_mutations_with_tolerance(list1, list2, tolerance=2):
    matched_mutations = []

    for mutation1 in list1:
        amino_acid1, position1 = mutation1.split()[1], int(mutation1.split()[-1])

        for mutation2 in list2:
            amino_acid2, position2 = mutation2.split()[1], int(mutation2.split()[-1])

            # Confronto sia dell'amminoacido che della posizione con tolleranza
            if amino_acid1 == amino_acid2 and abs(position1 - position2) <= tolerance:
                matched_mutations.append(mutation1)
                break

    return matched_mutations