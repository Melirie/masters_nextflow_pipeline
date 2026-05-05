params {

    raw_url: String
    raw_counts: Path
    donor_category: ArrayList
    donor_disease: ArrayList
    organ_groups: ArrayList
    cell_type: String
    batch: String
    cov_of_interest: String
    comparison_key: String

}

// Workflow block

// Include modules
include { load_data } from './modules/load_data.nf'
include { subsetting_1 } from './modules/subsetting_1.nf'
// include { initial_qc } from './modules/initial_QC.nf'
// include { metadata_correction } from './modules/metadata_correction.nf'
include { preprocessing } from './modules/preprocessing.nf'
include { normalization; pca; neighbors; umap } from './modules/standard_sc_workflow.nf'
include { qc_plot } from './modules/QC_plot.nf'
include { subsetting_2 } from './modules/subsetting_2.nf'
include { sample_qc } from './modules/sample_qc.nf'
include { scvi } from './modules/scvi.nf'
include { mrvi } from './modules/mrVI.nf'
// include { kbet; }
include { differential_abundance_calc; differential_abundance_plot } from './modules/differential_abundance.nf'
include { umap_after_batch_cor } from './modules/umap_after_batch_cor.nf'
include { umap_plot } from './modules/umap_plot.nf'
include { differential_expression_calc; differential_expression_plot } from './modules/differential_expression.nf'


workflow {
    main:
    load_data(params.raw_url, 'adata_all.h5ad')
    subsetting_1(load_data.out, params.donor_category, params.donor_disease, params.organ_groups, 'adata_subset.h5ad')
    // metadata_correction(subsetting_1.out, 'sex_metadata_correction_plot.pdf','adata_subset_metadata_fixed.h5ad')
    // initial_qc(metadata_correction.out.h5ad, 'adata_qc.h5ad')
    preprocessing(subsetting_1.out, 'sex_metadata_correction_plot.pdf','adata_subset_metadata_fixed.h5ad')
    qc_plot(preprocessing.out.h5ad, 'qc_plot.png')
    normalization(preprocessing.out.h5ad, 'adata_subset_normalized.h5ad')
    pca(normalization.out, 'adata_subset_pca.h5ad')
    neighbors(pca.out, 'adata_subset_neighbors.h5ad')
    umap(neighbors.out, 'umap_unintegrated.pdf','adata_clean.h5ad')
    subsetting_2(umap.out.h5ad, params.cell_type, "adata_${params.cell_type}.h5ad")
    sample_qc(subsetting_2.out, "adata_${params.cell_type}_clean.h5ad")
    // output_array = [subsetting_2.out, umap.out]
    // output_ch = channel.of(output_array)
    scvi(sample_qc.out, params.batch, "scvi_model_${params.cell_type}")
    mrvi(sample_qc.out, params.batch, params.cov_of_interest, "mrvi_model_${params.cell_type}")
    umap_after_batch_cor(mrvi.out, scvi.out, subsetting_2.out, "adata_all_batch_${params.cell_type}.h5ad")
    bases_ch = channel.of('unintegrated', 'scvi', 'mrvi_u', 'mrvi_z').view()
    plot_inputs = umap_after_batch_cor.out.combine(bases_ch)
    umap_plot(plot_inputs)
    differential_abundance_calc(mrvi.out, umap_after_batch_cor.out, params.comparison_key, "da_${params.cell_type}.nc")
    differential_abundance_plot(mrvi.out, umap_after_batch_cor.out, differential_abundance_calc.out, params.comparison_key, "da_${params.cell_type}.pdf")
    differential_expression_calc(mrvi.out, umap_after_batch_cor.out, params.comparison_key, "deg_${params.cell_type}.nc")

    publish:
    loaded_data = load_data.out
    subset_data = subsetting_1.out
    // metadata_data = metadata_correction.out.h5ad
    // sex_correction_plot = metadata_correction.out.plot
    qc_data = preprocessing.out.h5ad
    qc_plot = qc_plot.out
    umap_plot = umap.out.plot
    sex_correction_plot = preprocessing.out.plot
    adata_clean = umap.out.h5ad
    adata_cell_type_subset = subsetting_2.out
    scvi_model = scvi.out
    mrvi_model = mrvi.out
    mrvi_adata = umap_after_batch_cor.out
    umap_after_bc = umap_plot.out
    da_results = differential_abundance_calc.out
    da_plot = differential_abundance_plot.out
    deg_results = differential_expression_calc.out

}

output {

    loaded_data {
        path { "../data" }
    }

    subset_data {
        path { "../data" }
    }

    qc_data {
        path { "../data" }
    }

    qc_plot {
        path { "./plots" }
    }

    umap_plot {
        path { "./plots" }
    }

    sex_correction_plot {
        path { "./plots" }
    }

    adata_clean {
        path { "../data" }
    }

    adata_cell_type_subset {
        path { "../data" }
    }

    scvi_model {
        path { "./models" }
    }

    mrvi_model {
        path { "./models" }
    }

    mrvi_adata {
        path { "../data" }
    }

    umap_after_bc {
        path { "./plots" }
    }

    da_results {
        path { "./tables" }
    }

    da_plot {
        path { "./plots" }
    }

    deg_results {
        path { "/mnt/data/melina" }
    }


}