process sample_qc {
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

# 3. Identify samples that have 200 or more cells
keep_samples = sample_counts[sample_counts >= 200].index

# 4. Subset the adata to only include those samples
adata = adata[adata.obs['sampleID'].isin(keep_samples)].copy()

# Optional: Print how many samples were removed
removed = len(sample_counts) - len(keep_samples)
print(f"Removed {removed} samples with less than 200 cells.")
print(f"Remaining cells: {adata.n_obs}")

adata.write('${output_file}')

END_PYTHON
    """
}