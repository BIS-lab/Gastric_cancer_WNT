import pandas as pd
import seaborn as sns
import matplotlib as plt 
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib import rcParams
from scipy.stats import wilcoxon

### Visualization

# Load data
normalized_expr_data = pd.read_csv("./tmm_normalized_counts_all_sample_final.csv",index_col="Unnamed: 0")
metadata = pd.read_csv("./2_metadata.csv")
tx2gene = pd.read_csv("./0_result/star_salmon/salmon_tx2gene.tsv", sep="\t",header=None)
tx2gene.columns = ["Transcript_id","Gene_id","Gene_name"]



valid_gene_id = []
for i in ["ERBB2","WNT7B"]:
    gene_tmp = tx2gene.loc[tx2gene["Gene_name"] == i,"Gene_id"].iloc[0]
    valid_gene_id.append(gene_tmp)

snv_expr = normalized_expr_data.loc[valid_gene_id, :].copy()
snv_expr["gene_id"] = snv_expr.index

df_melted = snv_expr.melt(id_vars="gene_id", var_name="Source", value_name="Expression")
df_melted["Source"] = df_melted["Source"].apply(lambda x: metadata[metadata["sample"]==x]["Source"].values[0])
epsilon = 1
df_melted["Expression"] = np.log2(df_melted["Expression"] + epsilon)

# samples & colors
stripplot_samples = ["OO9", "OO14", "OO66", "DD191"]
boxplot_samples = ["gastric cancer", "normal stomach"]
legend_labels = stripplot_samples + boxplot_samples
set2_palette = sns.color_palette("Set2", len(legend_labels))
legend_handles = [
    mpatches.Patch(color=color, label=label)
    for color, label in zip(set2_palette, legend_labels)
]


fig, axs = plt.subplots(1, 2, figsize=(8, 5))
rcParams['pdf.fonttype'] = 42

for i, gene_id in enumerate(valid_gene_id):
    ax = axs[i]

    # data extraction
    df_temp = df_melted[df_melted["gene_id"] == gene_id].copy()
    df_temp["Source"] = pd.Categorical(
        df_temp["Source"],
        categories=legend_labels,
        ordered=True
    )

    # boxplot (bulk sample)
    df_box = df_temp[~df_temp["Source"].isin(stripplot_samples)]
    sns.boxplot(
        x="Source", y="Expression", data=df_box,
        palette="Set2", ax=ax, width=0.6
    )

    # stripplot (PDO sample)
    df_strip = df_temp[df_temp["Source"].isin(stripplot_samples)]
    sns.stripplot(
        x="Source", y="Expression", data=df_strip,
        palette="Set2", ax=ax, size=6, jitter=0.3,
        marker="o", linewidth=0.3
    )


    gene_name = tx2gene.loc[tx2gene["Gene_id"] == gene_id, "Gene_name"].iloc[0]
    ax.set_title(f"{gene_name}", fontsize=10)
    ax.set_ylabel("log2 TMM normalized Expression", fontsize=9)
    ax.tick_params(axis='x', rotation=45)
    for label in ax.get_xticklabels():
        label.set_ha('right')


plt.tight_layout()
plt.subplots_adjust(wspace=0.4)
plt.savefig("Expression_comparison_plot.pdf", bbox_inches='tight')





### Statistics

results = []

for gene in ["ERBB2", "WNT7B"]:
    
    gene_id = tx2gene.loc[tx2gene["Gene_name"] == gene, "Gene_id"].iloc[0]
    
    # Extract expression
    GC_group_data = df_melted.loc[(df_melted["gene_id"] == gene_id) & 
                                  (df_melted["Source"] == 'gastric cancer'), "Expression"]
    Normal_group_data = df_melted.loc[(df_melted["gene_id"] == gene_id) & 
                                      (df_melted["Source"] == 'normal stomach'), "Expression"]

    # Single sample Wilcoxon test
    for ssample in ["OO9", "OO14", "OO66", "DD191"]:
        single_sample = float(df_melted.loc[(df_melted["gene_id"] == gene_id) & 
                                            (df_melted["Source"] == ssample), "Expression"])
        
        # compare to GC
        differences_GC = [single_sample - x for x in GC_group_data]
        stat_GC, p_value_GC = wilcoxon(differences_GC, alternative="greater", zero_method="pratt") # greater, less, two-sided
        
        # Compare to Normal
        differences_Normal = [single_sample - x for x in Normal_group_data]
        stat_NM, p_value_NM = wilcoxon(differences_Normal, alternative="greater")
        print(stat_NM, p_value_NM)
        
        results.append([gene, ssample, "Gastric Cancer", stat_GC, p_value_GC])
        results.append([gene, ssample, "Normal Stomach", stat_NM, p_value_NM])

df_results = pd.DataFrame(results, columns=["Gene", "Sample", "Comparison Group", "Wilcoxon Statistic", "P-Value"])

df_results.to_csv("Wilcoxon_Test_Results_greater.csv", index=False)
