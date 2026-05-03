process metadata_correction {
    input:
    path input_file
    val plot_output_file
    val output_file

    output:
    path "${plot_output_file}", emit: plot
    path "${output_file}", emit: h5ad

    script:


    """
    python3 <<-END_PYTHON

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.api.types import CategoricalDtype

# 1. Load data
adata = sc.read_h5ad('${input_file}')

# 2. Fix Unknown Sex Categories & Generate Plot
sex_lower = adata.obs['sex'].astype(str).str.lower()
is_female = sex_lower == 'female'
is_male = sex_lower == 'male'
is_unknown = sex_lower.isin(['unknown', 'nan', 'none'])
genes_to_use = [g for g in ['XIST', 'RPS4Y1', 'EIF1AY', 'DDX3Y'] if g in adata.var_names]

def get_flat_mean(subset_adata, genes):
    return np.ravel(subset_adata[:, genes].X.mean(axis=0))

avg_female = get_flat_mean(adata[is_female], genes_to_use)
avg_male = get_flat_mean(adata[is_male], genes_to_use)

ref_df = pd.DataFrame([avg_female, avg_male], 
                      columns=genes_to_use, 
                      index=['REF: Average Female', 'REF: Average Male'])

if is_unknown.any():
    adata_unk = adata[is_unknown].copy()
    if hasattr(adata_unk.obs['donorID_unified'], 'cat'):
        adata_unk.obs['donorID_unified'] = adata_unk.obs['donorID_unified'].cat.remove_unused_categories()
    
    unk_df = adata_unk[:, genes_to_use].to_df()
    unk_df['donorID_unified'] = adata_unk.obs['donorID_unified'].values
    unk_means = unk_df.groupby('donorID_unified', observed=True).mean()
    final_df = pd.concat([ref_df, unk_means])
else:
    final_df = ref_df

# Reshape and plot
plot_df = final_df.reset_index().melt(id_vars='index', var_name='Gene', value_name='Mean_Expression')
plot_df.rename(columns={'index': 'Donor'}, inplace=True)

plt.figure(figsize=(18, 8))
sns.set_style("white")
sns.barplot(data=plot_df, x='Donor', y='Mean_Expression', hue='Gene', palette={'XIST': '#e63946', 'RPS4Y1': '#457b9d', 'EIF1AY': '#1d3557', 'DDX3Y': '#a8dadc'})
plt.title("Sex Marker Profiles", fontsize=16, fontweight='bold')
plt.axvline(x=1.5, color='black', linestyle='--', alpha=0.4)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("${plot_output_file}")


# 3. Age & Metadata corrections
# Match age of specific donors
adata.obs.loc[adata.obs["donorID_unified"].isin(["D143", "D144", "D145"]), "age_unified"] = "47-80"

# Change donorID_unified of D146
new_labels = ["D146A", "D146B"]
if not isinstance(adata.obs["donorID_unified"].dtype, pd.CategoricalDtype):
    adata.obs["donorID_unified"] = adata.obs["donorID_unified"].astype('category')
adata.obs["donorID_unified"] = adata.obs["donorID_unified"].cat.add_categories(new_labels)

adata.obs.loc[adata.obs["donor_id"] == "S279", "donorID_unified"] = "D146A"
adata.obs.loc[adata.obs["donor_id"] == "S280", "donorID_unified"] = "D146B"
adata.obs["donorID_unified"] = adata.obs["donorID_unified"].cat.remove_unused_categories()

# Assign "0-3" to donor D146B
age_order = ["0-3", "4-7", "9-12", "13-17", "18-34", "35-54", "47-80", "55-74", "75+"]
adata.obs["age_unified"] = adata.obs["age_unified"].astype(CategoricalDtype(categories=age_order, ordered=True))
adata.obs.loc[adata.obs["donorID_unified"] == "D146B", "age_unified"] = "0-3"

# 4. Generate donor_disease_category
adata.obs["donor_disease_category"] = (adata.obs["donor_disease"].astype(str) + "_" + adata.obs["donor_category"].astype(str))
adata.obs["donor_disease_category"] = adata.obs["donor_disease_category"].replace({"organ_donor_control": "healthy_control", "control_control": "healthy_control"})
adata.obs["donor_disease_category"] = adata.obs["donor_disease_category"].astype("category").cat.add_categories('pediatric_healthy_control')

pediatric_mask = (adata.obs['donor_disease_category'] == 'healthy_control') & (adata.obs['age_unified'].isin(['0-3', '4-7', '9-12', '13-17']))
adata.obs.loc[pediatric_mask, 'donor_disease_category'] = 'pediatric_healthy_control'

# 5. Rename gene indices with safety check
if not adata.var['gene_symbols'].is_unique:
    raise ValueError("Gene symbols are not unique! Cannot set as index.")

adata.var.set_index('gene_symbols', inplace=True)

# 6. Save
adata.write('${output_file}')

END_PYTHON
    """
}