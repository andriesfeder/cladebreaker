process ANALYZE_REPORT {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11':
        'quay.io/biocontainers/python:3.11' }"

    input:
    path(reports)

    output:
    path "*.analyze_report.txt" , emit: report
    path "*.analyze_report.tsv" , emit: results
    path "*.analyze_report.json", emit: json
    path "*_mqc.tsv"            , emit: mqc
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'cladebreaker'
    """
    analyze_report.py \\
        --alpha ${params.gsi_alpha} \\
        ${args} \\
        --out-prefix ${prefix} \\
        ${reports}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^Python //')
    END_VERSIONS
    """
}
