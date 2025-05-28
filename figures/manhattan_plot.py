import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import uniform

data = pd.read_csv("gwas_results.csv")

# Manhattan plot
def manhattan_plot(data, genome_wide_line=5e-8):
    data['-log10(P)'] = -np.log10(data['P'])
    data['ind'] = range(len(data))
    data_grouped = data.groupby('CHR')

    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ['#4C72B0', '#55A868']
    x_labels = []
    x_labels_pos = []

    for i, (chr_num, group) in enumerate(data_grouped):
        group.plot(kind='scatter', x='ind', y='-log10(P)', color=colors[i % len(colors)], ax=ax, s=10)
        x_labels.append(chr_num)
        x_labels_pos.append((group['ind'].iloc[-1] + group['ind'].iloc[0]) / 2)

    ax.axhline(-np.log10(genome_wide_line), color='red', linestyle='--', lw=1)
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(P-value)')
    ax.set_title('Manhattan Plot')
    sns.despine()
    plt.tight_layout()
    plt.savefig("manhattan_plot.png", dpi=300)
    plt.show()
