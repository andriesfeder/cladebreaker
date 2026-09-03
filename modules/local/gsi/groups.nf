process GSI_GROUPS {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::dendropy=5.0.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dendropy:5.0.13--pyhdfd78af_0':
        'quay.io/biocontainers/dendropy:5.0.13--pyhdfd78af_0' }"

    input:
    path(tree)
    path(query_ids)
    path(reference_ids)

    output:
    path "gsi_groups.tsv", emit: groups
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gsi_groups.py \\
        --tree ${tree} \\
        --query-ids ${query_ids} \\
        --reference-ids ${reference_ids} \\
        --query-label ${params.gsi_query_label} \\
        --reference-label ${params.gsi_reference_label} \\
        ${args} \\
        --out gsi_groups.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dendropy: \$(python -c "import dendropy; print(dendropy.__version__)")
    END_VERSIONS
    """
}
