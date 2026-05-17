process qc_plot {
    
    input:
    path input_file

    output:
    path "${input_file.baseName}_qc.png"

    script:
    """
    python3 <<-END_PYTHON

import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad('${input_file}')

# --- QC PLOTTING BY STUDY ---
# Standard Scanpy violin plot for QC metrics
with plt.rc_context({'figure.figsize': (12, 8)}):
    axes = sc.pl.violin(
        adata, 
        ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], 
        groupby='study', 
        multi_panel=True,
        stripplot=False,
        show=False
    )
    # Customize the layout if needed
    plt.savefig('${input_file.baseName}_qc.png', bbox_inches='tight')
    plt.close()

END_PYTHON
    """
}