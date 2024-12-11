import matplotlib.pyplot as plt
# Unique Alignments
# Definire i dati
labels = ['Rand', 'Rand Fr', 'Rand PAM30', 'RITA S', 'RITA M', 'RITA L', 'RITA XL',
          'SpikeGPT2', 'SARITA S', 'SARITA M', 'SARITA L', 'SARITA XL']
values = [0, 0, 0.05, 0, 0, 0, 0, 0.12, 17, 58.70, 72.92, 60.8]

# Creare un grafico a barre
plt.figure(figsize=(10, 6))
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
          '#7f7f7f', '#bcbd22', '#17becf', '#ff9896', '#98df8a']
# Creare un grafico a barre
plt.barh(labels, values, color=colors)

# Aggiungere etichette sui valori con un font più grande
for index, value in enumerate(values):
    plt.text(value + 0.2, index, f'{value:.2f}%', va='center', fontsize=12)

# Aggiungere etichette e titolo con un font più grande
plt.xlabel('Percentage', fontsize=14)
plt.title('Unique Alignment Performance Across Various Models', fontsize=16)

# Invertire l'asse Y per avere i valori dall'alto verso il basso
plt.gca().invert_yaxis()

# Removing borders
for spine in plt.gca().spines.values():
    spine.set_visible(False)

# Modificare la dimensione delle etichette degli assi
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Aggiungere una griglia
plt.grid(axis='x', linestyle='--', alpha=0.7)

# Aggiustare i limiti dell'asse X per dare risalto ai valori piccoli
plt.xlim(0, 80)

# Salvare il grafico come immagine
plt.tight_layout()
plt.savefig('Unique_Alignments')
plt.show()

# Definire i nuovi dati mantenendo lo stesso ordine e nomi esatti dell'asse Y come nel grafico precedente
labels_corrected = ['Rand', 'Rand Fr', 'Rand PAM30', 'RITA S', 'RITA M', 'RITA L', 'RITA XL',
                    'SpikeGPT2', 'SARITA S', 'SARITA M', 'SARITA L', 'SARITA XL']
#values_corrected = [0, 0, 10.1, 0, 0, 0, 0, 12.5, 98.8, 91.8, 89.2, 89.2]
values_corrected = [99.9, 99.9, 99.9, 99.9, 99.9, 99.9, 99.9, 12.5, 98.9, 99.8, 99.5, 99.7]

# Creare un nuovo grafico a barre orizzontale con i nomi corretti sulle colonne dell'asse Y
plt.figure(figsize=(10, 6))
plt.barh(labels_corrected, values_corrected, color=colors)

# Aggiungere etichette sui valori con un font più grande
for index, value in enumerate(values_corrected):
    plt.text(value + 0.2, index, f'{value:.1f}%', va='center', fontsize=12)

# Aggiungere etichette e titolo con un font più grande
plt.xlabel('Percentage', fontsize=14)
plt.title('Coverage of Test Set Mutations', fontsize=16)

# Invertire l'asse Y per mantenere l'ordine
plt.gca().invert_yaxis()

# Removing borders
for spine in plt.gca().spines.values():
    spine.set_visible(False)

# Modificare la dimensione delle etichette degli assi
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Aggiungere una griglia
plt.grid(axis='x', linestyle='--', alpha=0.7)

# Aggiustare i limiti dell'asse X
plt.xlim(0, 100)

# Visualizzare il grafico
plt.tight_layout()
plt.savefig('test_coverage')
plt.show()

# Updated colors for better distinction
updated_colors = [
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],  # Blu pastello per IN_TRAIN_AND_TEST, Rosa pastello per IN_TEST_ONLY, Verde pastello per IN_TRAIN_ONLY, Giallo pastello per NOVEL
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A']
    ]

# Data for the charts based on the provided pie charts
data_list = [
    {"labels": ["In training and test", "In test only", "In training only", "Novel"],
     "percentages": [0.39, 0.08, 0.05, 0.48],
     "colors": updated_colors[2]},

    {"labels": ["In training and test", "In test only", "In training only", "Novel"],
     "percentages": [0.30, 0.11, 0.04, 0.55],
     "colors": updated_colors[3]},

    {"labels": ["In training and test", "In test only", "In training only", "Novel"],
     "percentages": [0.31, 0.11, 0.04, 0.54],
     "colors": updated_colors[4]},

    {"labels": ["In training and test", "In test only", "In training only", "Novel"],
     "percentages": [0.31, 0.11, 0.04, 0.54],
     "colors": updated_colors[5]}
]

# Setting figure size larger for better visibility
fig, ax = plt.subplots(figsize=(14, 10))

# New order of indices for S, M, L, XL
new_order = [0, 1, 2, 3]

# Setting up a larger width for each bar
width = 0.5
x_positions = range(len(data_list))

# Plotting stacked bars for each chart in the new order with increased width
for idx, new_idx in enumerate(new_order):
    bottom = 0
    for i, p in enumerate(data_list[new_idx]['percentages']):
        ax.bar(idx, p, bottom=bottom, width=width, color=updated_colors[new_idx][i], edgecolor='black')
        bottom += p
        if p != 0:
            ax.text(idx, bottom - p / 2, f'{p:.2f}', ha='center', va='center', fontsize=16)

# Removing borders
for spine in plt.gca().spines.values():
    spine.set_visible(False)

# Setting x-axis labels in the new order
ax.set_xticks(range(len(new_order)))
ax.set_xticklabels(['SARITA-S', 'SARITA-M', 'SARITA-L', 'SARITA-XL'], fontsize=16)

# Removing borders
for spine in plt.gca().spines.values():
    spine.set_visible(False)


# Adjusting y-axis to a common scale
ax.set_ylim(0, 1)
ax.set_xlabel('Models', fontsize=16)
ax.set_ylabel('Percentage', fontsize=16)
ax.set_title('Comparison of SARITA Models in Mutation Generation', fontsize=18)
# Adjusting the size of the tick labels on both axes
ax.tick_params(axis='both', which='major', labelsize=14)
# Definisci i nuovi nomi per la legenda
custom_legend_labels = ['In training and test', 'In test only', 'In training only', 'Novel']

# Aggiungi la legenda personalizzata
ax.legend(custom_legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, prop={'size': 14})

plt.tight_layout()
plt.savefig('Comparison_Models.png')
plt.show()



import matplotlib.pyplot as plt

# Colori aggiornati per migliorare la distinzione
updated_colors = [
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],  # Blu pastello per IN_TRAIN_AND_TEST, Rosa pastello per IN_TEST_ONLY, Verde pastello per IN_TRAIN_ONLY, Giallo pastello per NOVEL
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A'],
    ['#A3C1DA', '#F7A1A1', '#A8E6A3', '#FDE68A']
]

# Dati per i grafici a barre
data_list = {
    "Rand": {"percentages": [0.30, 0.11, 0.04, 0.55], "colors": updated_colors[4]},
    "Rand Fr": {"percentages":[0.30, 0.11, 0.04, 0.55], "colors": updated_colors[4]},
    "Rand PAM30": {"percentages": [0.18, 0.08, 0.02, 0.72], "colors": updated_colors[4]},
    "RITA S": {"percentages": [0.30, 0.11, 0.04, 0.55], "colors": updated_colors[4]},
    "RITA M": {"percentages": [0.30, 0.11, 0.04, 0.55], "colors": updated_colors[4]},
    "RITA L": {"percentages": [0.30, 0.11, 0.04, 0.55], "colors": updated_colors[4]},
    "RITA XL": {"percentages": [0.31, 0.11, 0.04, 0.55], "colors": updated_colors[4]},  # Dati mancanti
    "SpikeGPT2": {"percentages": [1, 0, 0, 0], "colors": updated_colors[5]},  # Dati mancanti
    "SARITA S": {"percentages": [0.39, 0.05, 0.08, 0.48], "colors": updated_colors[2]},  # Dati mancanti
    "SARITA M": {"percentages": [0.30, 0.11, 0.04, 0.55], "colors": updated_colors[1]},  # Dati mancanti
    "SARITA L": {"percentages": [0.31, 0.11, 0.04, 0.54], "colors": updated_colors[0]},  # Dati mancanti
    "SARITA XL": {"percentages": [0.31, 0.11, 0.04, 0.54], "colors": updated_colors[3]}  # Dati mancanti
}

# Colore per i modelli mancanti
missing_colors = ['#d3d3d3','#FFFFFF']

# Definire l'ordine corretto dei modelli
labels_new_order = ['Rand', 'Rand Fr', 'Rand PAM30', 'RITA S', 'RITA M', 'RITA L', 'RITA XL',
                    'SpikeGPT2', 'SARITA S', 'SARITA M', 'SARITA L', 'SARITA XL']

# Impostazione delle dimensioni della figura
fig, ax = plt.subplots(figsize=(16, 10))

width = 0.5  # Larghezza delle barre

# Aggiunta delle barre impilate per ciascun grafico nel nuovo ordine
for idx, model in enumerate(labels_new_order):
    bottom = 0
    model_data = data_list.get(model)
    if model_data and model_data['percentages'] is not None:
        # Modelli con dati
        for i, p in enumerate(model_data['percentages']):
            ax.bar(idx, p, bottom=bottom, width=width, color=model_data['colors'][i], edgecolor='black')
            bottom += p
            if p != 0:
                ax.text(idx, bottom - p / 2, f'{p:.2f}', ha='center', va='center', fontsize=16)
    else:
        # Modelli mancanti in grigio con asterisco
        ax.bar(idx, 50, bottom=bottom, width=width, color=missing_colors[1], edgecolor='white')
        ax.text(idx, 0.5, '*', ha='center', va='center', fontsize=22, color='black')

# Impostazione delle etichette sull'asse X nel nuovo ordine
ax.set_xticks(range(len(labels_new_order)))
ax.set_xticklabels(labels_new_order, fontsize=16)

# Removing borders
for spine in plt.gca().spines.values():
    spine.set_visible(False)

# Regolazione dell'asse Y per una scala comune
ax.set_ylim(0, 1)
ax.set_xlabel('Models', fontsize=16)
ax.set_ylabel('Fraction', fontsize=16)
ax.set_title('Comparison of Models in Mutation Generation', fontsize=18)

# Regolazione della dimensione delle etichette degli assi
ax.tick_params(axis='both', which='major', labelsize=14)

# Aggiunta della legenda con il colore grigio chiaro
legend_labels = ['IN_TRAIN_AND_TEST', 'IN_TEST_ONLY', 'IN_TRAIN_ONLY', 'NOVEL']
handles = [plt.Rectangle((0,0),1,1, color=updated_colors[0][0]),
           plt.Rectangle((0,0),1,1, color=updated_colors[0][1]),  # Rosso per IN_TEST_ONLY
           plt.Rectangle((0,0),1,1, color=updated_colors[0][2]),  # Verde per IN_TRAIN_ONLY
           plt.Rectangle((0,0),1,1, color=updated_colors[0][3]),
           plt.Rectangle((0,0),1,1, color=missing_colors[0])]  # Grigio chiaro per sequenze non valide
ax.legend(handles, legend_labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=5, prop={'size': 14})

plt.tight_layout()

# Salvare il grafico
plt.savefig('Comparison_Models_Distribution.png')

# Visualizzare il grafico
plt.show()
