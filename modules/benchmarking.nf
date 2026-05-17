process benchmarking {
    input:
    path input_adata
    val output_plot

    output:
    path "*.svg", emit: plot_dir
    path "benchmarking_results.csv", emit: table

    script:


    """
    
    python3 <<-END_PYTHON

import numpy as np
import scanpy as sc
import pandas as pd

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

adata = sc.read_h5ad("${input_adata}")

import time

distance_keys = [
    "unintegrated_distances", 
    "neighbors_scvi_distances", 
    "neighbors_mrvi_u_distances", 
    "neighbors_mrvi_z_distances"
]

biocons = BioConservation(isolated_labels=True)

start = time.time()
bm = Benchmarker(
    adata,
    batch_key="donorID_unified",
    label_key="level_3_annot",
    embedding_obsm_keys=["X_pca", "X_scvi", "X_mrvi_u", "X_mrvi_z"],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=biocons,
    batch_correction_metrics=BatchCorrection(),
    n_jobs=-1,
)

bm.benchmark()
end = time.time()
print(f"Time: {int((end - start) / 60)} min {int((end - start) % 60)} sec")

bm.plot_results_table(save_dir = ".")
df = bm.get_results(min_max_scale=False)
df.to_csv("benchmarking_results.csv")

END_PYTHON
    """
}
