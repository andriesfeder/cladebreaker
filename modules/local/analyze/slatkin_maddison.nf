process SLATKIN_MADDISON {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::dendropy=5.0.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dendropy:5.0.13--pyhdfd78af_0':
        'quay.io/biocontainers/dendropy:5.0.13--pyhdfd78af_0' }"

    input:
    path(tree)
    path(groups)

    output:
    path "*.slatkin_maddison.tsv" , emit: results
    path "*.slatkin_maddison.json", emit: json
    path "*_mqc.tsv", emit: mqc
    
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: 'cladebreaker'
    def rooting   = params.gsi_root == 'outgroup' ? "--root outgroup --outgroup ${params.gsi_outgroup}" : "--root ${params.gsi_root}"
    def subset    = params.gsi_groups_to_test ? "--groups-to-test ${params.gsi_groups_to_test}" : ''
    def unlabeled = params.gsi_ignore_unlabeled ? '--ignore-unlabeled' : ''
    def perms     = "--permutations ${params.gsi_permutations} --seed ${params.gsi_seed}"
    """
    slatkin_maddison.py \\
        --tree ${tree} \\
        --groups ${groups} \\
        ${rooting} \\
        ${subset} \\
        ${unlabeled} \\
        ${perms} \\
        ${args} \\
        --out-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dendropy: \$(python -c "import dendropy; print(dendropy.__version__)")
    END_VERSIONS
    """
}
