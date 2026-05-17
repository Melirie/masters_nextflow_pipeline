process subsetting_1 {
    input:
    path input_file
    val donor_category
    val donor_disease
    val organ_groups
    val output_file

    output:
    path "${output_file}"

    script:

    // Convert Groovy lists to comma-separated strings
    def str_cat = donor_category.join(',')
    def str_dis = donor_disease.join(',')
    def str_org = organ_groups.join(',')

    """
    
    python3 -c "

import scanpy as sc

adata = sc.read_h5ad('${input_file}')

categories = '${str_cat}'.split(',')
diseases   = '${str_dis}'.split(',')
organs     = '${str_org}'.split(',')

# Clean up steps

# Clear Pairwise Matrices
if hasattr(adata, 'obsp'):
    for key in list(adata.obsp.keys()):
        del adata.obsp[key]

# Clear Metadata
for key in ['neighbors', 'umap', 'leiden']:
    if key in adata.uns:
        del adata.uns[key]

# Clear old embeddings
for key in ['X_scANVI', 'X_umap']:
    if key in adata.obsm:
        del adata.obsm[key]

# Clear old uns that are not needed anymore
stale_keys = ['neighbors', 'umap', 'leiden', 'pca', 'log1p']
for key in stale_keys:
    if key in adata.uns:
        del adata.uns[key]

#adata = adata.copy() # Breaking the view, had issues here before but maybe check if still necessary

# Subsetting for groups defined in the params
# Select conditions: Diseases, Donor Category (for excluding fetal tissue), Organ Groups
keep_cells = (adata.obs['donor_disease'].isin(diseases)) & \
             (adata.obs['donor_category'].isin(categories)) & \
             (adata.obs['organ_groups'].isin(organs))

#adata = adata[keep_cells].copy()
adata = adata[keep_cells]

adata.write_h5ad('${output_file}', compression=None)
"
    """
}