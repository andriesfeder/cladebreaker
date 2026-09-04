process ANALYZE_PDF {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::matplotlib-base=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/matplotlib:3.1.2--2':
        'quay.io/biocontainers/matplotlib:3.1.2--2' }"

    input:
    path(reports)
    path(tree)
    path(groups)

    output:
    path "*.analyze_report.pdf", emit: pdf
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'cladebreaker'
    """
    analyze_pdf.py \\
        --tree ${tree} \\
        --groups ${groups} \\
        ${args} \\
        --out-prefix ${prefix} \\
        ${reports}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
