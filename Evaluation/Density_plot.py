import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Estrai i dati dal CSV

data = pd.read_csv('/Users/utente/Desktop/Varcovid/GenSeq/Model/Evaluation/calc/merged_substitution_info.csv')

# Creazione del dizionario subs_dict
subs_dict = {
    'train': data[data['group'] == 'Train']['substitution'],
    'test': data[data['group'] == 'Test']['substitution'],
    'predict': data[data['group'] == 'Predict']['substitution']
}

# Creazione della lista di dizionari df_dicts con gestione degli errori di conversione
df_dicts = []
for k, v in subs_dict.items():
    for x in v:
        try:
            # Filtrare e convertire a numero intero
            numeric_value = int(''.join(filter(str.isdigit, x)))
            df_dicts.append({'group': k, 'value': numeric_value})
        except ValueError:
            # Gestione di valori non convertibili
            print(f"Valore non numerico trovato e ignorato: {x}")

# Creazione del DataFrame plot_df
plot_df = pd.DataFrame(df_dicts)

REFSEQ_START = 0
palette = sns.color_palette('colorblind', len(plot_df['group'].unique()))

# Creiamo una figura e assi
plt.figure(figsize=(10, 6))

# Tracciamo il grafico KDE per ciascun gruppo senza riempimento
for group in plot_df['group'].unique():
    subset = plot_df[plot_df['group'] == group]
    sns.kdeplot(subset.loc[subset['value'] >= REFSEQ_START + 14, 'value'],
                label=group, bw_adjust=0.2)

# Aggiungiamo etichette e titolo
plt.xlabel('Value')
plt.ylabel('Density')
plt.title('Kernel Density Estimate (KDE) Plot')
plt.legend(title='Group')
plt.show()

