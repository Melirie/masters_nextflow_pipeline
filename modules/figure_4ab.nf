process figure_4ab {
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
# Fixed: Ensured this is a list to avoid the 'KeyError: l' string iteration error
plot_columns = ['study', 'donor_disease_category']

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
        
        # Handle case where subplots might not be a 2D array if ROWS=1, COLS=1
        if PLOTS_PER_PAGE > 1:
            axes_flat = axes.flatten()
        else:
            axes_flat = [axes]
        
        # Determine columns for this specific page
        start_idx = p * PLOTS_PER_PAGE
        page_columns = plot_columns[start_idx : start_idx + PLOTS_PER_PAGE]
        
        for i, col in enumerate(page_columns):
            ax = axes_flat[i]
            
            # Use Scanpy to plot onto the specific axis
            # legend_loc=None prevents Scanpy from drawing a duplicate legend
            sc.pl.embedding(
                adata, 
                basis=current_basis, 
                color=col, 
                ax=ax, 
                show=False,
                legend_loc=None,
                title=f"{col}_${base}"
            )
            
            # Position the legend to the right of the axis (outside the plot)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(
                handles, 
                labels,
                loc='upper left',
                bbox_to_anchor=(1.02, 1),
                fontsize='5', 
                frameon=False,
                markerscale=0.5
            )
            
            # Clean up plot aesthetics
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.title.set_fontsize(8)

        # Hide any empty subplots on the final page
        for j in range(len(page_columns), PLOTS_PER_PAGE):
            axes_flat[j].axis('off')

        # Use rect to ensure tight_layout makes room for the legends on the right
        plt.tight_layout(rect=[0, 0, 0.9, 1])
        
        pdf.savefig(fig)
        plt.close(fig)

print(f"Successfully created {output_filename}")

END_PYTHON
    """
}