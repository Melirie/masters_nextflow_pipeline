process differential_abundance_calc {
    input:
    path input_model
    path input_adata
    val comparison_key
    val output_file

    output:
    path "${output_file}"

    script:


    """
    
    python3 <<-END_PYTHON

import scanpy as sc
import scvi
import torch
import os
import matplotlib.pyplot as plt
import time
from scvi.external import MRVI
import pandas as pd
import numpy as np
import xarray as xr



# 1. Configuration
# We pull the path from the Bash variable
input_file = os.environ.get('DATA_PATH')

# 1. Load your "Fresh" data
adata = sc.read_h5ad("${input_adata}")

# 2. Load Model C
model_mrvi = MRVI.load("${input_model}", adata=adata)
adata

print(f"Starting differential abundance calcultation for {adata.n_obs} cells...")


de_results = model_mrvi.differential_abundance(
    adata=adata,                       
    sample_cov_keys=["${comparison_key}"], 
    batch_size=512,
    compute_log_enrichment=True
)

print(f' calculation complete.')


# 1. Cast categorical coordinates/data to strings to avoid TypeError
for var in de_results.coords:
    if isinstance(de_results[var].dtype, pd.CategoricalDtype):
        de_results[var] = de_results[var].astype(str)

for var in de_results.data_vars:
    if isinstance(de_results[var].dtype, pd.CategoricalDtype):
        de_results[var] = de_results[var].astype(str)


# 3. Save the xarray Dataset
de_results.to_netcdf("${output_file}")
print(f"Results successfully saved")


END_PYTHON
    """
}


process differential_abundance_plot {
    input:
    path input_model
    path input_adata
    path input_da
    val comparison_key
    val output_file

    output:
    path "${output_file}"

    script:


    """
    
    python3 <<-END_PYTHON

import scanpy as sc
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PathCollection

# 1. Load Data
adata = sc.read_h5ad("${input_adata}")
da_res = xr.open_dataset("${input_da}")

# 2. Define Baseline Logic
# We map target categories to their specific baselines
control_map = {
    "CD_disease": "CD_control",
    "UC_disease": "UC_control",
    "PIBD_disease": "pediatric_healthy_control"
}
# Global baseline
global_ref = "healthy_control"

# Get all unique categories in the data
categories = da_res["${comparison_key}"].values

# 3. Initialize PDF
with PdfPages('${output_file}') as pdf:
    
    # --- Page 1: Reference Annotations ---
    sc.pl.embedding(
        adata, basis='umap_mrvi_u', color=['level_2_annot', 'level_3_annot'], 
        frameon=False, ncols=2, show=False 
    )
    for ax in plt.gcf().axes:
        for artist in ax.collections:
            if isinstance(artist, PathCollection): artist.set_rasterized(True)
    pdf.savefig(bbox_inches='tight')
    plt.close()

    # --- Loop through categories and apply your logic ---
    for cat in categories:
        # Determine which baselines to compare against
        baselines = []
        
        # 1. Everything gets compared to healthy_control (if cat isn't healthy_control itself)
        if cat != global_ref:
            baselines.append(global_ref)
        
        # 2. Respective disease comparison (if mapped in control_map)
        specific_ref = control_map.get(cat)
        if specific_ref and specific_ref in categories:
            baselines.append(specific_ref)

        for ref in baselines:
            if cat == ref: continue
            
            print(f"Plotting: {cat} vs {ref}")
            
            # Calculate Log Prob Ratio
            target_probs = da_res.donor_disease_category_log_probs.loc[{"${comparison_key}": cat}]
            ref_probs = da_res.donor_disease_category_log_probs.loc[{"${comparison_key}": ref}]
            lfc = target_probs - ref_probs
            
            # Map to adata
            adata.obs["DA_lfc"] = lfc.to_series().reindex(adata.obs_names).values

            # Plot
            sc.pl.embedding(
                adata, basis='umap_mrvi_u', color=['DA_lfc'], 
                cmap='coolwarm', vmin=-1, vmax=1, size=2.0, frameon=False,
                title=f"LFC: {cat} vs {ref}", show=False
            )
            
            # Manual rasterization for performance
            for ax in plt.gcf().axes:
                for artist in ax.collections:
                    if isinstance(artist, PathCollection): artist.set_rasterized(True)
            
            pdf.savefig(bbox_inches='tight', dpi=300)
            plt.close()

END_PYTHON
    """
}