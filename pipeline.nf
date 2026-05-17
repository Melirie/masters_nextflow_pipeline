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
include { preprocessing } from './modules/preprocessing.nf'
include { standard_sc_workflow } from './modules/standard_sc_workflow.nf'
include { qc_plot } from './modules/QC_plot.nf'
include { subsetting_2 } from './modules/subsetting_2.nf'
include { sample_qc } from './modules/sample_qc.nf'
include { scvi } from './modules/scvi.nf'
include { mrvi } from './modules/mrVI.nf'
include { differential_abundance_calc; differential_abundance_plot } from './modules/differential_abundance.nf'
include { umap_after_batch_cor } from './modules/umap_after_batch_cor.nf'
// include { umap_plot } from './modules/umap_plot.nf'
include { differential_expression_calc; differential_expression_plot } from './modules/differential_expression.nf'
include { benchmarking } from './modules/benchmarking.nf'
include { figure_2_plots } from './modules/figure_2.nf'
include { figure_3a } from './modules/figure_3a.nf'



workflow {
    main:
    load_data(params.raw_url, 'adata_all.h5ad')
    subsetting_1(load_data.out, params.donor_category, params.donor_disease, params.organ_groups, 'adata_subset.h5ad')
    preprocessing(subsetting_1.out, 'sex_metadata_correction_plot.pdf','adata_subset_metadata_fixed.h5ad')
    standard_sc_workflow(preprocessing.out.h5ad,'adata_clean.h5ad','umap_unintegrated.pdf')
    subsetting_2(standard_sc_workflow.out.h5ad, params.cell_type, "adata_${params.cell_type}.h5ad")
    // sample_qc(subsetting_2.out, "adata_${params.cell_type}_clean.h5ad")
    // output_array = [subsetting_2.out, umap.out]
    // output_ch = channel.of(output_array)
    scvi(subsetting_2.out, params.batch, "scvi_model_${params.cell_type}")
    mrvi(subsetting_2.out, params.batch, params.cov_of_interest, "mrvi_model_${params.cell_type}")
    umap_after_batch_cor(mrvi.out, scvi.out, subsetting_2.out, "adata_all_batch_${params.cell_type}.h5ad")
    // files_to_qc = ["/mnt/data/melina/adata_clean.h5ad", "/mnt/data/melina/adata_all_batch_Epithelial.h5ad"]
    // This creates a channel that emits each file as a Path object
    // ch_qc = Channel.fromPath(files_to_qc)
    // Parallel execution starts here
    qc_plot(umap_after_batch_cor.out)
    bases_ch = channel.of('unintegrated', 'scvi', 'mrvi_u', 'mrvi_z')
    plot_inputs = umap_after_batch_cor.out.combine(bases_ch)
    figure_3a(plot_inputs)
    differential_abundance_calc(mrvi.out, umap_after_batch_cor.out, params.comparison_key, "da_${params.cell_type}.nc")
    differential_abundance_plot(mrvi.out, umap_after_batch_cor.out, differential_abundance_calc.out, params.comparison_key, "da_${params.cell_type}.pdf")
    differential_expression_calc(mrvi.out, umap_after_batch_cor.out, params.comparison_key, "deg_${params.cell_type}_uc_control_baseline.nc")
    // differential_expression_plot(differential_expression_calc.out, umap_after_batch_cor.out, "deg_${params.cell_type}.pdf")
    benchmarking(umap_after_batch_cor.out, "benchmarking_${params.cell_type}")
    figure_2_plots(load_data.out, subsetting_1.out, subsetting_2.out)

    publish:
    loaded_data = load_data.out
    qc_plot = qc_plot.out
    umap_plot = standard_sc_workflow.out.plot
    sex_correction_plot = preprocessing.out.plot
    adata_umap = standard_sc_workflow.out.h5ad
    adata_cell_type_subset = subsetting_2.out
    scvi_model = scvi.out
    mrvi_model = mrvi.out
    mrvi_adata = umap_after_batch_cor.out
    umap_after_bc = figure_3a.out
    da_results = differential_abundance_calc.out
    da_plot = differential_abundance_plot.out
    deg_results = differential_expression_calc.out
    // deg_plot = differential_expression_plot.out
    benchmarking_plot = benchmarking.out.plot_dir
    benchmarking_table = benchmarking.out.table
    sex_correction_heatmap = figure_2_plots.out.sex_heatmap
    subset_table = figure_2_plots.out.table
    age_plot = figure_2_plots.out.age_plot
    figure_3 = figure_3a.out

}

output {

    loaded_data {
        path { "/mnt/data/melina" }
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

    adata_umap {
        path { "/mnt/data/melina" }
    }

    adata_cell_type_subset {
        path { "/mnt/data/melina" } // after subsetting 2
    }

    scvi_model {
        path { "/mnt/data/melina/models" }
    }

    mrvi_model {
        path { "/mnt/data/melina/models" }
    }

    mrvi_adata {
        path { "/mnt/data/melina" }
    }

    umap_after_bc {
        path { "./plots" }
    }

    da_results {
        path { "/mnt/data/melina" }
    }

    da_plot {
        path { "./plots" }
    }

    deg_results {
        path { "/mnt/data/melina" }
    }

    benchmarking_plot {
        path { "./plots" }
    }

    benchmarking_table {
        path { "./tables" }
    }

    sex_correction_heatmap {
        path { "./plots" }
    }

    subset_table {
        path { "./tables" }
    }

    age_plot {
        path { "./plots" }
    }

    figure_3 {
        path { "./plots" }
    }


}