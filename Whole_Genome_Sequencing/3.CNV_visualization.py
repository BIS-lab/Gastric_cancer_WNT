import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams

wnt_genes = ['WNT1', 'WNT2', 'WNT2B', 'WNT3', 'WNT4', 'WNT5A', 'WNT5B', 'WNT6', 'WNT7A', 'WNT7B', 'WNT8A', 'WNT8B', 'WNT9A', 'WNT9B', 'WNT10A', 'WNT10B', 'WNT11', 'WNT16',
             'TP53', 'ERBB2', 'ERBB3', 'FGFR2', 'CDH1', 'GATA6', 'KRAS', 'RNF43', 'MYC']

data_dir = "./CNV_visualizaiton"

all_data = pd.DataFrame()

rcParams['pdf.fonttype'] = 42


for file_name in os.listdir(data_dir):
    if file_name.endswith("gene.txt"):
        file_path = os.path.join(data_dir, file_name)
        sample_name = file_name.split(".")[0].split("_gene")[0]  
        df = pd.read_csv(file_path, sep='\t')  
        all_data[sample_name] = np.nan
        for target_gene in wnt_genes:
            df_gene = df[df['gene'] == target_gene]  
            if not df_gene.empty:  
                all_data.loc[target_gene, sample_name] = 2*2**(float(df_gene["log2"].values[0]))

all_data = all_data.astype(float)
custom_sample_order = ['GA353T',
 'GA372T',
 'YPGC-057', 
 'YPGC-075', 
 'YPGC-082', 
 'YPGC-090', 
 'YPGC-105', 
 'YPGC-162', 
 'OO9',
 'OO14',
 'OO66',
 'DD191']
all_data = all_data.reindex(columns=custom_sample_order)


plt.figure(figsize=(10, 20))
heatmap = sns.heatmap(all_data, cmap="coolwarm", annot=True, center=2, vmax=4, fmt=".2f", cbar_kws={'label': 'Copy Number'})


for text in heatmap.texts:
    value = float(text.get_text()) 
    if value >= 3.5:
        text.set_color('yellow') 
        text.set_weight('bold')
    elif value >= 2*2**0.2:
        text.set_color('red')  
        text.set_weight('bold')
    elif value <= 2*2**-0.2:
        text.set_color('blue') 
plt.title(f"Copy Number Variation")


plt.savefig("Figure_S6E.pdf", format='pdf', bbox_inches='tight')  # PDF로 저장
plt.show()

plt.savefig("Figure_S6E.png", format='png', dpi=300, bbox_inches='tight')  # PNG로 저장 (dpi는 해상도)
plt.show()