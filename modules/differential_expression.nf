process differential_expression_calc {
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



# Hardware Setup
# torch.set_float32_matmul_precision('high')


# 1. Configuration
# We pull the path from the Bash variable
input_file = os.environ.get('DATA_PATH')
job_id = os.environ.get('SLURM_JOB_ID', 'standalone')

# 1. Load your "Fresh" data
adata = sc.read_h5ad("${input_adata}")

# 2. Load Model C
model_mrvi = MRVI.load("${input_model}", adata=adata) #final mrvi model from job 171, sample_key = sample_ID, batch_key = donor
adata



current_categories = model_mrvi.sample_info["donor_disease_category"].cat.categories.tolist()
print(current_categories)
print("Reordering covariate keys...")


sample_cov_keys = adata.obs["donor_disease_category"]
model_mrvi.sample_info["donor_disease_category"] = model_mrvi.sample_info["donor_disease_category"].cat.reorder_categories(
    ['healthy_control', 'pediatric_healthy_control', 'PIBD_disease', 'UC_control', 'UC_disease']

)

print(f"Starting deg calcultation for {adata.n_obs} cells...")

de_results = model_mrvi.differential_expression(
    adata=adata,
    sample_subset=None,                       
    sample_cov_keys=["donor_disease_category"],
    batch_size=4,                     # Default is 
    use_vmap=False,                    
    mc_samples=5,                     # monte carlo samples was 15
    filter_inadmissible_samples=True,  # filter out-of-distribution samples prior to performing the analysis
    add_batch_specific_offsets=False,    # Whether to offset the design matrix by adding batch-specific offsets to the design matrix
    store_lfc=True                        # Whether to store log fold changes in the results
)


# Save the xarray Dataset
de_results.to_netcdf("${output_file}", engine="netcdf4")



END_PYTHON
    """
}


process differential_expression_plot {
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


END_PYTHON
    """
}