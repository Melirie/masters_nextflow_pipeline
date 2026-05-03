process umap_plot {
    input:
    path input_file
    val basis
    val output_plot

    output:
    path "${output_plot}"

    script:
    """
    python3 <<-END_PYTHON

import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

adata = sc.read_h5ad('${input_file}')

# Define the columns you want to plot
plot_columns = ['level_1_annot', 'study', 'donor_disease_category', 'organ_groups']

# Open a multi-page PDF file
with PdfPages('${output_plot}') as pdf:
    for col in plot_columns:
        # Create the plot for this specific column
        sc.pl.embedding(
            adata, 
            basis="umap_${basis}_u", 
            color=col, 
            show=False
        )
        
        # Save the current figure to the PDF
        pdf.savefig(bbox_inches='tight')
        
        # Crucial: Close the figure so the next one starts fresh
        plt.close()

END_PYTHON
    """
}