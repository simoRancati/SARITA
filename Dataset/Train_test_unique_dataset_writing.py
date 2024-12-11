import statistics as st
from utils import *

Continenti = ['World']
for l in Continenti:
    print('starting reading files')
    # Read the File csv
    print('Read file csv')
    #/ Users / utente / Desktop / Varcovid / DeepAutoCov / raw_data / metadata_tsv_2023_11_06 / metadata.csv
    df = pd.read_csv("/blue/salemi/share/varcovid/set_2024_08_19/metadata.csv") #
    print('End Read File csv')

    print('Start Read the FASTA files')
    sequences, headers = read_protein_sequences_header("/blue/salemi/share/varcovid/set_2024_08_19/spikeprot0818/spikes.fasta")
    print('End Read FASTA File')

    # Merging Data
    print('Mearging Data')
    df.iloc[:, 4] = df.iloc[:, 4].astype(type(headers[0]))
    df_ordinato = pd.DataFrame(headers, columns=['ID']).merge(df, left_on='ID', right_on=df.columns[4])
    df_ordinato['Sequences'] = sequences
    df = df_ordinato.drop('ID', axis=1)
    sequences = df["Sequences"].values.tolist()
    df = df.drop(['Sequences'], axis=1)
    metadata_nation = df.to_numpy()
    print('End Mearging Data')

    # In the Fasta format, the asterisk is used to indicate the end of a sequence. The ending asterisk is removed in some sequences.
    print("\033[1mRemoval of the final asterisk in the Fasta format.\033[0m")
    sequences_nation = [remove_asterisks(s) for s in sequences]

    Dimension = len(sequences_nation)
    print('\033[1mNumber of spike protein is: ' + str(Dimension) + '\033[0m')
    # Lunghezza delle sequenze
    Length = []
    for sequence in sequences_nation:
        Length.append(len(sequence))

    print('\033[1mFilter sequences with length < 1000 \033[0m')
    sequences_filtered_min_1000 = [x for i, x in enumerate(sequences_nation) if Length[i] >= 1000]
    index_1 = [i for i, x in enumerate(Length) if x >= 1000]
    print('\033[1mAggiorno il file Metadata\033[0m')
    metadata_filtered_min_1000 = metadata_nation[index_1]
    Dimension_fil_min_1000 = len(sequences_filtered_min_1000)
    print('The number of spike proteins after eliminating sequences with length less than 1000 is: ' + str(Dimension_fil_min_1000))

    print('\033[1mCalculation of filtered lengths less than 1000\033[0m')
    Length_filtered_min_1000 = []
    for sequence in sequences_filtered_min_1000:
        Length_filtered_min_1000.append(len(sequence))

    # Seleziono le Sequenze Valide
    print('\033[1mEvaluate the valid sequences contained in the database\033[0m')
    valid_sequences, invalid_sequences, valid_indices, invalid_indices = validate_sequences(sequences_filtered_min_1000)
    print('\033[1mEvaluate : \033[0m')
    print('There are '+str(len(valid_sequences))+' valid sequences in database')
    print('There are '+str(len(invalid_sequences))+' invalid sequences in the database')

    print('\033[1mupdate filtered \033[0m')
    metadata_valid_indices = metadata_filtered_min_1000[valid_indices]
    metadata_valid_invalid_indices = metadata_filtered_min_1000[invalid_indices]

    print('\033[1mCalculate the new length of the filtered and valid sequences \033[0m')
    Length_filtered_min_1000_valid = []
    for sequence in valid_sequences:
        Length_filtered_min_1000_valid.append(len(sequence))

    print(st.median(Length_filtered_min_1000_valid))

    print('\033[1mFilter sequences by the length within the median\033[0m')
    extreme_up = 1269 #st.median(Length_filtered_min_1000_valid) - 30
    extreme_low = 1275 #st.median(Length_filtered_min_1000_valid) + 30
    index_1,sequences_valid = filter_sequences(valid_sequences, extreme_up, extreme_low)
    metadata_valid_indices_length = metadata_valid_indices[index_1]
    print('Selected ' + str(len(sequences_valid)) + ' with length between ' + str(extreme_up) + ' and ' + str(extreme_low))

    print("\033[1mFilter sequences by dates \033[0m")
    metadata_off,sequences_off, metadata_not_off, sequences_not_off = filter_row_by_column_length_sostitution(metadata_valid_indices_length, sequences_valid, 5, 10) #Valuto che le date siano valide
    # In NLP applications the data should not be repeated so many
    # times because the data would be biased for this reason I am going
    # to select only a few
    index_seq = select_indices(sequences_off)
    sequences_off = [sequences_off[i] for i in index_seq]
    metadata_off = metadata_off[index_seq, :]

    print("\033[1mThe number of sequences filtered with dates is :\033[0m"+str(len(metadata_off)))

    print("\033[1mReorder metadata file \033[0m")
    metadata_tot = insert_sequence_as_column(metadata_off,metadata_off[:,5],sequences_off)

    sequences = list(metadata_tot[:,24])
    metadata = metadata_tot[:, 0:23]

    print(' Weeks')
    indices_by_week = split_weeks_custom(metadata[:,5],2024,8,18)
    print(len(indices_by_week))
    seq_week=[]
    for i in range(0,len(indices_by_week)):
        seq_week.append(len(indices_by_week[i]))
    print(seq_week)

    write_csv_dataset(metadata, l)

    for i in range(0,len(indices_by_week)):
        indices = indices_by_week[i]
        sequences_for_week = []
        identifier_for_week = []
        week = i+1
        # Creating Dataset
        os.makedirs('/blue/salemi/share/varcovid/GenSeq/dataset_Aug_2024_fasta_'+l + '/' + str(week))
        for m,index in enumerate(indices):
            sequences_for_week.append(sequences[index])
            identifier_for_week.append(metadata[index,4])
            create_multiple_fasta_file_fast(identifier_for_week, sequences_for_week, 'sequences.fasta', week, l)

