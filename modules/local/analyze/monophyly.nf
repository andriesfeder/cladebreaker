process MONOPHYLY {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::dendropy=5.0.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dendropy:5.0.13--pyhdfd78af_0':
        'quay.io/biocontainers/dendropy:5.0.13--pyhdfd78af_0' }"

    input:
    path(tree)
    path(groups)

    output:
    path "*.monophyly.tsv" , emit: results
    path "*.monophyly.json", emit: json
    path "*_mqc.tsv", emit: mqc
    path "*_clade_breaking.tsv", emit: clade_breaking
    path "*.rooted.nwk"        , emit: rooted_tree
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: 'cladebreaker'
    def rooting   = params.gsi_root == 'outgroup' ? "--root outgroup --outgroup ${params.gsi_outgroup}" : "--root ${params.gsi_root}"
    def subset    = params.gsi_groups_to_test ? "--groups-to-test ${params.gsi_groups_to_test}" : ''
    def unlabeled = params.gsi_ignore_unlabeled ? '--ignore-unlabeled' : ''

    """
    monophyly.py \\
        --tree ${tree} \\
        --groups ${groups} \\
        ${rooting} \\
        ${subset} \\
        ${unlabeled} \\
        ${args} \\
        --out-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dendropy: \$(python -c "import dendropy; print(dendropy.__version__)")
    END_VERSIONS
    """
}
