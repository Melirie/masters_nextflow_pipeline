process figure_2_plots {
    input:
    path input_file_1
    path input_file_2
    path input_file_3

    output:
    path "tables_cell_counts.csv", emit: table
    path "age_distribution_cells.svg", emit: age_plot
    path "sex_correction_heatmap.svg", emit: sex_heatmap

    script:
    """
    python3 <<-END_PYTHON
    import scanpy as sc
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    # --- PART 1: Summary Table ---
    file_paths = ["${input_file_1}", "${input_file_2}", "${input_file_3}"]
    dataset_names = ["Initial data", "Subsetting 1", "Subsetting 2"] 

    results = []
    for path, name in zip(file_paths, dataset_names):
        adata = sc.read_h5ad(path, backed='r')
        stats = {
            'Dataset': name,
            'Total Cells': adata.n_obs,
            'Total Unique Genes': adata.n_vars,
            'Total Samples': adata.obs['sampleID'].nunique() if 'sampleID' in adata.obs.columns else "N/A",
            'Total Organ Groups': adata.obs['organ_groups'].nunique() if 'organ_groups' in adata.obs.columns else "N/A",
            'Total Donors': adata.obs['donorID_unified'].nunique() if 'donorID_unified' in adata.obs.columns else "N/A"
        }
        results.append(stats)

    summary_df = pd.DataFrame(results)
    summary_t = summary_df.set_index('Dataset').T
    summary_t.to_csv("tables_cell_counts.csv", index=True)

    # --- PART 2: Age Distribution (Using Input 2) ---
    adata = sc.read_h5ad("${input_file_2}")
    
    if 'age_unified' in adata.obs.columns:
        age_counts = adata.obs['age_unified'].value_counts(dropna=False).sort_index()
        age_counts.index = [str(i) if pd.notna(i) else "Missing/NA" for i in age_counts.index]
        colors = ['#e63946' if x == "Missing/NA" else 'skyblue' for x in age_counts.index]

        plt.figure(figsize=(12, 6))
        age_counts.plot(kind='bar', color=colors, edgecolor='black', width=0.8)
        plt.title("Number of Cells per Age Group (including NA)", fontsize=14, fontweight='bold')
        plt.xlabel("Age Group", fontsize=12)
        plt.ylabel("Cell Count", fontsize=12)
        
        for i, v in enumerate(age_counts.values):
            plt.text(i, v + (v * 0.01), f"{v:,}", ha='center', fontweight='bold')

        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig("age_distribution_cells.svg")
        plt.close()
    else:
        # Create dummy plot if column missing to avoid Nextflow output error
        plt.figure().savefig("age_distribution_cells.svg")

    # --- PART 3: Heatmap (Using Input 2) ---
    if 'gene_symbols' in adata.var.columns:
        adata.var_names = adata.var['gene_symbols'].astype(str)
        adata.var_names_make_unique()

    genes = ['XIST', 'RPS4Y1']
    available_genes = [g for g in genes if g in adata.var_names]
    
    if len(available_genes) > 0 and 'donorID_unified' in adata.obs.columns:
        matrix = adata[:, available_genes].X.toarray() if hasattr(adata.X, 'toarray') else adata[:, available_genes].X
        final_df = pd.DataFrame(matrix, index=adata.obs['donorID_unified'], columns=available_genes)
        heatmap_df = final_df.groupby(level=0).mean()

        plt.figure(figsize=(14, 6))
        sns.heatmap(heatmap_df.T, cmap="RdBu_r", center=0, linewidths=.5, cbar_kws={'label': 'Mean Expression'})
        plt.title("Sex Marker Expression Profile", fontsize=16, fontweight='bold')
        plt.xlabel("Donor", fontsize=12)
        plt.ylabel("Gene", fontsize=12)
        plt.tight_layout()
        plt.savefig("sex_correction_heatmap.svg")
        plt.close()
    else:
        plt.figure().savefig("sex_correction_heatmap.svg")

    END_PYTHON
    """
}