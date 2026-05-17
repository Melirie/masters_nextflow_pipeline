// This process includes normalisation, pca, neighbor calculation, umap calculation and plotting
process standard_sc_workflow {
    input:
    path input_file
    val output_file
    val output_plot

    output:
    path "${output_file}", emit: h5ad
    path "${output_plot}", emit: plot


    script:


    """
    
    python3 <<-END_PYTHON

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read_h5ad('${input_file}')

# Normalize library size to 10,000 reads per cell, create the "data" layer first by copying the counts
adata.layers["data"] = adata.layers["counts"].copy()

# Normalize the layer in-place
sc.pp.normalize_total(adata, target_sum=1e4, layer="data")

# Log transformation
sc.pp.log1p(adata, layer="data")
adata.layers["data"] = adata.X.copy()

# Check all layers
for layer_name, layer_data in adata.layers.items():
    print(f"Layer '{layer_name}': {layer_data.shape}")

# Check raw, which often causes a 'hidden' shape mismatch
if adata.raw is not None:
    print(f"Raw Object: {adata.raw.shape}")

# Calculate PCA
sc.tl.pca(adata)

# Calculate neighbors
sc.pp.neighbors(adata, use_rep="X_pca", key_added="unintegrated")

# Calculate UMAP
sc.tl.umap(adata, neighbors_key="unintegrated", key_added = "X_umap_unintegrated")
# adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"].copy()

# Write adata file into nf output
adata.write('${output_file}')

# Plot UMAP
sc.pl.embedding(
    adata, 
    basis='umap_unintegrated', 
    color='donor_disease_category',
    show=False
)

# Save plot into nf output
plt.savefig('${output_plot}', bbox_inches='tight')

END_PYTHON
    """
}