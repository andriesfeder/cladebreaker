process SNP_SEPARATION {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::dendropy=5.0.13" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dendropy:5.0.13--pyhdfd78af_0':
        'quay.io/biocontainers/dendropy:5.0.13--pyhdfd78af_0' }"

    input:
    path(distances)
    path(groups)

    output:
    path "*.snp_separation.tsv"      , emit: results
    path "*.snp_separation.json"     , emit: json
    path "*.snp_separation_pairs.tsv", emit: pairs
    path "*_mqc.tsv"                 , emit: mqc
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: 'cladebreaker'
    def subset    = params.gsi_groups_to_test ? "--groups-to-test ${params.gsi_groups_to_test}" : ''
    def unlabeled = params.gsi_ignore_unlabeled ? '--ignore-unlabeled' : ''
    """
    snp_separation.py \\
        --distances ${distances} \\
        --groups ${groups} \\
        --permutations ${params.gsi_permutations} \\
        --seed ${params.gsi_seed} \\
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
