process initial_qc {
    input:
    path input_file
    val output_file

    output:
    path "${output_file}"

    script:
    """
    python3 <<-END_PYTHON

import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('${input_file}')

# Filter samples with < 3 cells
counts = adata.obs['sampleID'].value_counts()
keep = counts[counts > 3].index
removed = counts[counts <= 3]

# Apply filter
adata = adata[adata.obs['sampleID'].isin(keep)].copy()

# 1. Save raw counts (from .raw.X) to a dedicated layer
adata.layers["counts"] = adata.raw.X.copy()

# Summary check
print(f"Raw counts layer now available: {list(adata.layers.keys())}")

n_genes_start = adata.n_vars
print(f"Starting with {n_genes_start} genes in active workspace.")

# We pull raw integers from .raw
adata.X = adata.raw[:, adata.var_names].X.copy()
print("Successfully swapped Raw Integers into adata.X.")

# Filter on raw counts
sc.pp.filter_genes(adata, min_counts=3)
n_genes_after_filter = adata.n_vars
genes_removed = n_genes_start - n_genes_after_filter

# We need to replace our raw slot, because the mismatch in dimensions would cause the batch corr to crash
adata.raw = adata
print("Successfully swapped raw.X and var with filtered counts")

print(f"Genes removed (min_counts=3 on RAW): {genes_removed}")

# Save raw counts also in counts layer as some tools expect this format
adata.layers["counts"] = adata.X.copy()


# compute QC statistics, we do not filter on them yet

# 1. Identify mitochondrial and ribosomal genes
# Standard Human prefixes: MT- for mitochondria, RPS/RPL for ribosome
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

# 2. Calculate QC metrics for both
# We use inplace=True to add 'pct_counts_mt' and 'pct_counts_ribo' to adata.obs
sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars=['mt', 'ribo'], 
    percent_top=None, 
    log1p=False, 
    inplace=True
)

# Load cell cycle genes from standard scanpy documentation
url = "https://raw.githubusercontent.com/scverse/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"
# We read it into a list, stripping any whitespace
raw_genes = pd.read_csv(url, header=None)[0].str.strip().tolist()

# check if present in dataset
s_genes = raw_genes[:43]
g2m_genes = raw_genes[43:]
cell_cycle_genes = [*s_genes, *g2m_genes]


s_genes_inter = [g for g in s_genes if g in adata.var_names]
g2m_genes_inter = [g for g in g2m_genes if g in adata.var_names]



# Use the intersecting genes we just verified
sc.tl.score_genes_cell_cycle(
    adata, 
    s_genes=s_genes_inter, 
    g2m_genes=g2m_genes_inter, 
    use_raw=False
)
# Quick check to ensure 'phase' exists now

adata.write("${output_file}")

END_PYTHON
    """
}