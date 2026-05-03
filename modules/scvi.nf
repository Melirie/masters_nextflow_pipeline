process scvi {
    input:
    path input_file
    val batch
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
import pandas as pd

# Hardware Setup
torch.set_float32_matmul_precision('high')

adata = sc.read_h5ad("${input_file}")

config = {
    "batch_size": 1024,
    "max_epochs": 1000,
    "early_stopping": True,
    "batch_keys": "${batch}"
}


# 3. Setup using the config
scvi.model.SCVI.setup_anndata(
    adata, 
    layer='counts', 
    batch_key=config["batch_keys"]
)

model = scvi.model.SCVI(adata) 


# 4. Train using the config
model.train(
    max_epochs=config["max_epochs"],
    batch_size=config["batch_size"],
    early_stopping=config["early_stopping"],
    enable_progress_bar=True, #testing if it works, can be changed to True for better monitoring
    accelerator='gpu', 
    devices=1
    # datasplitter_kwargs={
        # 'num_workers': 4,
        # 'pin_memory': True,
        # 'persistent_workers': True
    # }
)



# Save Model
model.save("${output_file}", overwrite=True)

END_PYTHON
    """
}