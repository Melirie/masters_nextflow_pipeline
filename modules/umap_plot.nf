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
import math

# 1. Load data
adata = sc.read_h5ad('${adata}')

# 2. Configuration
plot_columns = ['level_1_annot', 'level_2_annot', 'level_3_annot', 'study', 
                'donor_disease_category', 'organ_groups', 'sample_retrieval', 'assay']

current_basis = "umap_${base}"
output_filename = "umap_${base}_${params.cell_type}.pdf"

# A4 Dimensions and Grid Setup
A4_WIDTH = 8.27
A4_HEIGHT = 11.69
ROWS_PER_PAGE = 3
COLS_PER_PAGE = 2
PLOTS_PER_PAGE = ROWS_PER_PAGE * COLS_PER_PAGE

# 3. Plotting Loop
with PdfPages(output_filename) as pdf:
    num_pages = math.ceil(len(plot_columns) / PLOTS_PER_PAGE)
    
    for p in range(num_pages):
        # Initialize A4 Figure
        fig, axes = plt.subplots(ROWS_PER_PAGE, COLS_PER_PAGE, figsize=(A4_WIDTH, A4_HEIGHT))
        axes_flat = axes.flatten()
        
        # Determine columns for this specific page
        start_idx = p * PLOTS_PER_PAGE
        page_columns = plot_columns[start_idx : start_idx + PLOTS_PER_PAGE]
        
        for i, col in enumerate(page_columns):
            ax = axes_flat[i]
            
            # Use Scanpy to plot onto the specific axis
            sc.pl.embedding(
                adata, 
                basis=current_basis, 
                color=col, 
                ax=ax, 
                show=False,
                title=f"{col}_${base}"
            )
            
            # Shrink and position the legend
            # bbox_to_anchor puts the legend to the right of the plot
            leg = ax.legend(
                loc='best',
                bbox_to_anchor=(0.5, 0., 0.5, 0.5),
                fontsize='xx-small', 
                frameon=False,
                markerscale=0.4
            )
            
            # Clean up plot aesthetics
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.title.set_fontsize(8)

        # Hide any empty subplots on the final page
        for j in range(len(page_columns), PLOTS_PER_PAGE):
            axes_flat[j].axis('off')

        # Adjust layout to prevent overlap
        plt.tight_layout()
        
        pdf.savefig(fig)
        plt.close(fig)

print(f"Successfully created {output_filename}")

END_PYTHON
    """
}