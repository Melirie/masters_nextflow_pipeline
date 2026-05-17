process preprocessing {
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

adata = sc.read_h5ad("${input_file}")

# Raw counts handling
adata = adata.raw.to_adata().copy()
adata.layers["counts"] = adata.X.copy()

# Filter genes that are expressed in less than 10 cells
sc.pp.filter_genes(adata, min_counts=10)

# Define Sex genes
SEX_GENES = ['XIST', 'RPS4Y1', 'EIF1AY', 'DDX3Y']

for key in ['neighbors', 'umap', 'leiden', 'pca', 'log1p']:
    adata.uns.pop(key, None)
for key in ['X_scANVI', 'X_umap']:
    adata.obsm.pop(key, None)
adata.obsp.clear()
adata = adata.copy()

# 4. Metadata Fixes (Age & Donor IDs)
adata.obs.loc[adata.obs["donorID_unified"].isin(["D143", "D144", "D145"]), "age_unified"] = "47-80"

# Fix D146 pooled sample
if not isinstance(adata.obs["donorID_unified"].dtype, CategoricalDtype):
    adata.obs["donorID_unified"] = adata.obs["donorID_unified"].astype('category')
adata.obs["donorID_unified"] = adata.obs["donorID_unified"].cat.add_categories(["D146A", "D146B"])
adata.obs.loc[adata.obs["donor_id"] == "S279", "donorID_unified"] = "D146A"
adata.obs.loc[adata.obs["donor_id"] == "S280", "donorID_unified"] = "D146B"
adata.obs["donorID_unified"] = adata.obs["donorID_unified"].cat.remove_unused_categories()

# Age ordering
age_order = ["0-3", "4-7", "9-12", "13-17", "18-34", "35-54", "47-80", "55-74", "75+"]
adata.obs["age_unified"] = adata.obs["age_unified"].astype(CategoricalDtype(categories=age_order, ordered=True))
adata.obs.loc[adata.obs["donorID_unified"] == "D146B", "age_unified"] = "0-3"

# 5. Disease Category
# Create the combined column as a string first
adata.obs["donor_disease_category"] = (
    adata.obs["donor_disease"].astype(str) + 
    "_" + 
    adata.obs["donor_category"].astype(str)
)
print("Check old disease categories:")
print(adata.obs['donor_disease_category'].value_counts())

# Normalize the values first (This prepares the 'healthy_control' base)
rename_dict = {
    "organ_donor_control": "healthy_control",
    "control_control": "healthy_control"
}
adata.obs["donor_disease_category"] = adata.obs["donor_disease_category"].replace(rename_dict)

print("Check old disease categories:")
print(adata.obs['donor_disease_category'].value_counts())

# 3. convert to category and add your special label
adata.obs["donor_disease_category"] = adata.obs["donor_disease_category"].astype("category")
adata.obs["donor_disease_category"] = adata.obs["donor_disease_category"].cat.add_categories('pediatric_healthy_control')

# Finally, apply the pediatric filter
pediatric_ages = ['0-3', '4-7', '9-12', '13-17']
mask = (
    (adata.obs['donor_disease_category'] == 'healthy_control') & 
    (adata.obs['age_unified'].isin(pediatric_ages))
)

adata.obs.loc[mask, 'donor_disease_category'] = 'pediatric_healthy_control'


# 8. Gene Indexing
adata.var.index = adata.var['gene_symbols']

# 9. QC Metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], inplace=True)

# 10. Resolve Sex
new_female = ["D32", "D33", "D42", "D46", "D47", "D48", "D144"]
new_male = ["D31","D39", "D40", "D41", "D43", "D44", "D45", "D49", "D143", "D145"]
adata.obs.loc[adata.obs['donorID_unified'].isin(new_female), 'sex'] = 'female'
adata.obs.loc[adata.obs['donorID_unified'].isin(new_male), 'sex'] = 'male'

# --- DATA PREP FOR PLOT ---
# Calculate mean expression of markers per donor
sex_genes = ['XIST', 'RPS4Y1', 'EIF1AY', 'DDX3Y']
# Subset for genes that exist in data
valid_genes = [g for g in sex_genes if g in adata.var_names]

# Create temporary dataframe for calculation
temp_df = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X, 
    columns=adata.var_names, 
    index=adata.obs_names
)
temp_df['Donor'] = adata.obs['donorID_unified'].values
final_df = temp_df.groupby('Donor')[valid_genes].mean()

# --- PLOTTING ---
plot_df = final_df.reset_index().melt(id_vars='Donor', var_name='Gene', value_name='Mean_Expression')

plt.figure(figsize=(18, 8))
sns.set_style("white")
sns.barplot(
    data=plot_df, 
    x='Donor', 
    y='Mean_Expression', 
    hue='Gene', 
    palette={'XIST': '#e63946', 'RPS4Y1': '#457b9d', 'EIF1AY': '#1d3557', 'DDX3Y': '#a8dadc'}
)
plt.title("Sex Marker Profiles", fontsize=16, fontweight='bold')
plt.axvline(x=1.5, color='black', linestyle='--', alpha=0.4)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("${plot_output_file}")
plt.close()

# 6. Save
adata.write('${output_file}')

END_PYTHON
    """
}