process umap_plot {
    input:
    tuple path(adata), val(base)

    output:
    path "umap_${base}_${params.cell_type}.pdf"

    script:
    """
    python3 <<-END_PYTHON

import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

adata = sc.read_h5ad('${adata}')

# Define the columns you want to plot
plot_columns = ['level_1_annot', 'level_2_annot', 'level_3_annot', 'study', 'donor_disease_category', 'organ_groups', 'sample_retrieval', 'assay']

current_basis = "umap_${base}"
output_filename = "umap_${base}_${params.cell_type}.pdf"

# Open a multi-page PDF file
with PdfPages(output_filename) as pdf:
        for col in plot_columns:
            # basis looks for current_basis (e.g., 'umap_scvi') in adata.obsm
            sc.pl.embedding(
                adata, 
                basis=current_basis, 
                color=col, 
                title=f"{col} (${base})",
                show=False
            )
            
            pdf.savefig(bbox_inches='tight')
            plt.close()

END_PYTHON
    """
}