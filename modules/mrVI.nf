process mrvi {
    input:
    path input_file
    val batch
    val cov_of_interest
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

# Hardware Setup
torch.set_float32_matmul_precision('high')

config = {
    "batch_size": 1024,
    "max_epochs": 1000,
    "early_stopping": True,
    "sample_key": "${cov_of_interest}",
    "batch_key": "${batch}",
}

# 2. Data Loading
print("Loading data...")
adata = sc.read_h5ad("${input_file}")

sample_key = config["sample_key"]  # target covariate
batch_key=config["batch_key"]  # nuisance variable identifier
#labels_key=config["labels_key"]  # optional, for visualization
MRVI.setup_anndata(adata, layer="counts", sample_key=sample_key, batch_key=batch_key, backend="torch") #labels_key=labels_key
model = MRVI(adata)
 

# 4. Train using the config
model.train(
    max_epochs=config["max_epochs"],
    batch_size=config["batch_size"],
    early_stopping=config["early_stopping"],
    enable_progress_bar=True, #testing if it works, can be changed to True for better monitoring
    accelerator='gpu', 
    devices=1,
    datasplitter_kwargs={
        'num_workers': 8,
        'pin_memory': True,
        'persistent_workers': True
    }
)

model.save("${output_file}", overwrite=True)

END_PYTHON
    """
}