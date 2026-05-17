process figure_3a {
    input:
    tuple path(adata), val(base)

    output:
    path "fig3a_${base}_${params.cell_type}.pdf"

    script:
    """
    python3 <<-END_PYTHON

import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math

# 1. Load data
adata = sc.read_h5ad('${adata}')

# 2. Configuration
plot_columns = ['level_3_annot']

current_basis = "umap_${base}"
# Match exactly with Nextflow output: fig3a_...
output_filename = "fig3a_${base}_${params.cell_type}.pdf"

# A4 Dimensions and Grid Setup
A4_WIDTH = 8.27
A4_HEIGHT = 11.69
ROWS_PER_PAGE = 1
COLS_PER_PAGE = 1
PLOTS_PER_PAGE = ROWS_PER_PAGE * COLS_PER_PAGE

# 3. Plotting Loop
with PdfPages(output_filename) as pdf:
    num_pages = math.ceil(len(plot_columns) / PLOTS_PER_PAGE)
    
    for p in range(num_pages):
        fig, axes = plt.subplots(ROWS_PER_PAGE, COLS_PER_PAGE, figsize=(A4_WIDTH, A4_HEIGHT))
        
        if PLOTS_PER_PAGE > 1:
            axes_flat = axes.flatten()
        else:
            axes_flat = [axes]
        
        start_idx = p * PLOTS_PER_PAGE
        page_columns = plot_columns[start_idx : start_idx + PLOTS_PER_PAGE]
        
        for i, col in enumerate(page_columns):
            ax = axes_flat[i]
            
            # legend_loc='right margin' is the standard Scanpy way to show legends
            sc.pl.embedding(
                adata, 
                basis=current_basis, 
                color=col, 
                ax=ax, 
                show=False,
                legend_loc='right margin', 
                legend_fontsize='xx-small',
                title=f"{col}_${base}",
                frameon=False
            )
            
            # Clean up plot aesthetics
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.title.set_fontsize(8)

        # Hide any empty subplots on the final page
        for j in range(len(page_columns), PLOTS_PER_PAGE):
            axes_flat[j].axis('off')

        # Simple tight_layout with no rect constraints
        plt.tight_layout()
        
        pdf.savefig(fig)
        plt.close(fig)

print(f"Successfully created {output_filename}")

END_PYTHON
    """
}