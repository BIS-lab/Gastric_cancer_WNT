import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import spearmanr

# Load metadata and expression data
metadata = pd.read_csv("/vol/RZRZK/data/CCLE_metadata.csv")
expression_data = pd.read_csv("/vol/RZRZK/data/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

#Clean gene names
expression_data.columns = [col.split(" ")[0] if "(" in col else col for col in expression_data.columns]
expression_data.set_index('Unnamed: 0', inplace=True)
expression_data.index.name = "Cell line"
expression_data_transposed = expression_data.T

# Load gene list of interest
genelist = pd.read_csv("/vol/RZRZK/data/hallmark_kras_siganling.csv")

# Select the leftmost column (use the first column)
genelist = genelist["gene_symbol"].tolist()  # Remove duplicates and NaNs

wnt_gene = "WNT7B"

# Filter only valid genes
filtered_genes = [gene for gene in genelist if gene in expression_data_transposed.index]
if wnt_gene in filtered_genes:
    filtered_genes.remove(wnt_gene)  # Exclude WNT7B since it's the reference gene

# Remove non-cancer lines and clean data
data_transposed = expression_data_transposed.loc[:, expression_data_transposed.columns.isin(metadata["depmap_id"])]

# Calculate Pearson R and P-value
output_dir = "/vol/RZRZK/output/correlation"
os.makedirs(output_dir, exist_ok=True)

for type_of_column in ["Tissue_Type", "Cancer_Type", "Tissue", "Driver_type"]:
    all_tissue_correlation_results = {}
    all_tissue_pvalue_results = {}

    # Calculate correlation and p-value for all cell lines
    correlation_results = {}
    pvalue_results = {}

    for gene in filtered_genes:
        r_value, p_value = spearmanr(data_transposed.loc[wnt_gene], data_transposed.loc[gene])
        correlation_results[gene] = r_value if p_value < 0.05 else np.nan   # Keep only if p-value < 0.05
        pvalue_results[gene] = p_value

    all_tissue_correlation_results["Whole_Cell_lines"] = correlation_results
    all_tissue_pvalue_results["Whole_Cell_lines"] = pvalue_results

    # Calculate correlation by specific type
    for tissue in metadata[type_of_column].unique():
        tissue_specific_cell_lines = metadata.loc[metadata[type_of_column] == tissue, "depmap_id"]

        if len(tissue_specific_cell_lines) > 1:
            tissue_specific_data = data_transposed[tissue_specific_cell_lines]
            tissue_correlation_results = {}
            tissue_pvalue_results = {}

            for gene in filtered_genes:
                if gene in tissue_specific_data.index:
                    r_value, p_value = spearmanr(tissue_specific_data.loc[wnt_gene], tissue_specific_data.loc[gene])
                    tissue_correlation_results[gene] = r_value if p_value < 0.05 else np.nan  # Keep only if p-value < 0.05
                    tissue_pvalue_results[gene] = p_value
                else:
                    tissue_correlation_results[gene] = np.nan
                    tissue_pvalue_results[gene] = np.nan

            all_tissue_correlation_results[tissue] = tissue_correlation_results
            all_tissue_pvalue_results[tissue] = tissue_pvalue_results

    # Create DataFrames for correlation and p-values
    correlation_df = pd.DataFrame(all_tissue_correlation_results)
    pvalue_df = pd.DataFrame(all_tissue_pvalue_results)

    # Save CSV
    correlation_df.to_csv(f"{output_dir}/{type_of_column}_WNT7B_corr_hallmark_spearman.csv")
    
print(f"Filtered Spearman correlation results saved in '{output_dir}' directory")