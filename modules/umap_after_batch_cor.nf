process umap_after_batch_cor {
    input:
    path input_model
    path input_adata
    val output_file

    output:
    path "${output_file}"

    script:
    """
    python3 <<-END_PYTHON

import scanpy as sc
from scvi.external import MRVI

# 1. Load Data and Model
adata = sc.read_h5ad("${input_adata}")
model = MRVI.load("${input_model}", adata=adata)

# 2. Process the single model
label = "mrvi"
print(f"Processing Model {label}...")

# --- Latent Extraction ---
print("Getting latent representations...")
# u = sample-level (uncorrected), z = latent (corrected)
adata.obsm[f"X_{label}_u"] = model.get_latent_representation()
adata.obsm[f"X_{label}_z"] = model.get_latent_representation(give_z=True)

# --- Neighbors & UMAP for 'u' ---
print("Computing neighbors and UMAP for u...")
sc.pp.neighbors(adata, use_rep=f"X_{label}_u", key_added=f"neighbors_{label}_u")
sc.tl.umap(adata, neighbors_key=f"neighbors_{label}_u")
adata.obsm[f"X_umap_{label}_u"] = adata.obsm["X_umap"].copy()

# --- Neighbors & UMAP for 'z' ---
print("Computing neighbors and UMAP for z...")
sc.pp.neighbors(adata, use_rep=f"X_{label}_z", key_added=f"neighbors_{label}_z")
sc.tl.umap(adata, neighbors_key=f"neighbors_{label}_z")
adata.obsm[f"X_umap_{label}_z"] = adata.obsm["X_umap"].copy()

# Cleanup default UMAP slot to avoid downstream mixups
if "X_umap" in adata.obsm:
    del adata.obsm["X_umap"]

# 3. Save Output
adata.write_h5ad("${output_file}")


END_PYTHON
    """
}