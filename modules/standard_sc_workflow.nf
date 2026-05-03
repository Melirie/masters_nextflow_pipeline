process normalization {
    input:
    path input_file
    val output_file

    output:
    path "${output_file}"

    script:


    """
    
    python3 <<-END_PYTHON

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

adata = sc.read_h5ad('${input_file}')

adata.layers["counts"] = adata.X.copy()

# Normalize X for dimensionality reduction tools
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["data"] = adata.X.copy()
print("X and data are now normalized slots")

# Final check
print(f"Main Object (X): {adata.shape}")

# Check all layers
for layer_name, layer_data in adata.layers.items():
    print(f"Layer '{layer_name}': {layer_data.shape}")

# Check raw, which often causes the 'hidden' shape mismatch
if adata.raw is not None:
    print(f"Raw Object: {adata.raw.shape}")

adata.write('${output_file}')

END_PYTHON
    """
}


process pca {
    input:
    path input_file
    val output_file

    output:
    path "${output_file}"

    script:
    """
    python3 <<-END_PYTHON
import scanpy as sc


adata = sc.read_h5ad('${input_file}')

sc.tl.pca(adata)

adata.write('${output_file}')

END_PYTHON
    """
}

process neighbors {
    input:
    path input_file
    val output_file

    output:
    path "${output_file}"

    script:
    """
python3 <<-END_PYTHON
import scanpy as sc


adata = sc.read_h5ad('${input_file}')

sc.pp.neighbors(adata, use_rep="X_pca", key_added="unintegrated")

adata.write('${output_file}')

END_PYTHON
    """
}


process umap {
    input:
    path input_file
    val output_plot
    val output_file

    output:
    path "${output_plot}", emit: plot
    path "${output_file}", emit: h5ad

    script:
    """
    python3 <<-END_PYTHON

import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad('${input_file}')

sc.tl.umap(adata, neighbors_key="unintegrated")
adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"].copy()
del adata.obsm["X_umap"]

adata.write('${output_file}')

# 3. Create the Plot
sc.pl.embedding(
    adata, 
    basis='umap_unintegrated', 
    color='donor_disease_category',
    show=False
)

# 4. Save manually to the Nextflow output path
plt.savefig('${output_plot}', bbox_inches='tight')

END_PYTHON
    """
}