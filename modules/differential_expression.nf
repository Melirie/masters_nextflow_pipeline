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

scvi.settings.seed = 0

# 1. Load your "Fresh" data
adata = sc.read_h5ad("${input_adata}")

# 2. Load Model
model_mrvi = MRVI.load("${input_model}", adata=adata) #final mrvi model from job 171, sample_key = sample_ID, batch_key = donor
adata

#current_categories = model_mrvi.sample_info["donor_disease_category"].cat.categories.tolist()

sample_cov_keys = adata.obs["donor_disease_category"]
model_mrvi.sample_info["donor_disease_category"] = model_mrvi.sample_info["donor_disease_category"].cat.reorder_categories(
    #['healthy_control', 'pediatric_healthy_control', 'PIBD_disease', 'UC_control', 'UC_disease'] #if hc is supposed to be baseline
    ['UC_control', 'pediatric_healthy_control', 'PIBD_disease', 'healthy_control', 'UC_disease'] # if uc control is supposed to be baseline
)

de_results = model_mrvi.differential_expression(
    adata=adata,
    sample_subset=None,                       
    sample_cov_keys=["donor_disease_category"],
    batch_size=4,                     # below default
    use_vmap=False,                    
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