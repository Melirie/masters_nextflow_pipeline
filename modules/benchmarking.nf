process kbet {
    input:
    path input_model
    path input_adata
    val batch_key
    val output_file

    output:
    path "${output_file}"

    script:


    """
    
    python3 <<-END_PYTHON

import scib_metrics

batch_vector = adata.obs['${batch_key}'].values
results = scib_metrics.kbet(X, batches=batch_vector, alpha=0.05)



END_PYTHON
    """
}


process silhouette_batch {
    input:
    path input_model
    path input_adata
    val batch_key
    val output_file

    output:
    path "${output_file}"

    script:


    """
    
    python3 <<-END_PYTHON

import scib_metrics

batch_vector = adata.obs['${batch_key}'].values
results = scib_metrics.kbet(X, batches=batch_vector, alpha=0.05)



END_PYTHON
    """
}



process iLISI {
    input:
    path input_model
    path input_adata
    val batch_key
    val output_file

    output:
    path "${output_file}"

    script:


    """
    
    python3 <<-END_PYTHON

import scib_metrics

batch_vector = adata.obs['${batch_key}'].values
results = scib_metrics.kbet(X, batches=batch_vector, alpha=0.05)



END_PYTHON
    """
}


process PCR_comparison {
    input:
    path input_model
    path input_adata
    val batch_key
    val output_file

    output:
    path "${output_file}"

    script:


    """
    
    python3 <<-END_PYTHON

import scib_metrics

batch_vector = adata.obs['${batch_key}'].values
results = scib_metrics.kbet(X, batches=batch_vector, alpha=0.05)



END_PYTHON
    """
}