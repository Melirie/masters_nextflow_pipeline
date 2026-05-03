process subsetting_2 {
    input:
    path input_file
    val cell_type
    val output_file

    output:
    path "${output_file}"

    script:


    """
    
    python3 <<-END_PYTHON

import scanpy as sc

adata = sc.read_h5ad('${input_file}')

adata = adata[adata.obs["level_1_annot"]=='${cell_type}'].copy()
adata = adata[adata.obs["organ_groups"]=='Large_intestine'].copy()


sample_counts = adata.obs['sampleID'].value_counts()

# qc
# Identify samples that have 200 or more cells
keep_samples = sample_counts[sample_counts >= 200].index

# Subset the adata to only include those samples
adata = adata[adata.obs['sampleID'].isin(keep_samples)].copy()

# Print how many samples were removed
removed = len(sample_counts) - len(keep_samples)
print(f"Removed {removed} samples with less than 200 cells.")
print(f"Remaining cells: {adata.n_obs}")

adata.write_h5ad('${output_file}', compression=None)
END_PYTHON
    """
}