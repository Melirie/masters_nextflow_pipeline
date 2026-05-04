process umap_after_batch_cor {
    input:
    path input_model_mrvi
    path input_model_scvi
    path input_adata
    val output_file

    output:
    path "${output_file}"

    script:
    """
    python3 <<-END_PYTHON

import scanpy as sc
import scvi
from scvi.external import MRVI

# 1. Load Data and Model
adata = sc.read_h5ad("${input_adata}")

# Unintegrated PCA, neighbors, UMAP

sc.tl.pca(adata)
sc.pp.neighbors(adata, use_rep="X_pca", key_added="unintegrated")
sc.tl.umap(adata, neighbors_key="unintegrated")
adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"].copy()
del adata.obsm["X_umap"]

model_scvi = scvi.model.SCVI.load("${input_model_scvi}", adata=adata)
adata.obsm["X_scvi"] = model_scvi.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scvi", key_added="neighbors_scvi")
sc.tl.umap(adata, neighbors_key="neighbors_scvi")
adata.obsm["X_umap_scvi"] = adata.obsm["X_umap"].copy()
del adata.obsm["X_umap"]

model_mrvi = MRVI.load("${input_model_mrvi}", adata=adata)
adata.obsm[f"X_mrvi_u"] = model_mrvi.get_latent_representation()
adata.obsm[f"X_mrvi_z"] = model_mrvi.get_latent_representation(give_z=True)
sc.pp.neighbors(adata, use_rep=f"X_mrvi_u", key_added=f"neighbors_mrvi_u")
sc.tl.umap(adata, neighbors_key=f"neighbors_mrvi_u")
adata.obsm[f"X_umap_mrvi_u"] = adata.obsm["X_umap"].copy()
del adata.obsm["X_umap"]
sc.pp.neighbors(adata, use_rep=f"X_mrvi_z", key_added=f"neighbors_mrvi_z")
sc.tl.umap(adata, neighbors_key=f"neighbors_mrvi_z")
adata.obsm[f"X_umap_mrvi_z"] = adata.obsm["X_umap"].copy()
del adata.obsm["X_umap"]

# 3. Save Output
adata.write_h5ad("${output_file}")

END_PYTHON
    """
}