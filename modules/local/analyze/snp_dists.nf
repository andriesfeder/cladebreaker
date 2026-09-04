process SNP_DISTS {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::snp-dists=1.2.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snp-dists:1.2.0--h577a1d6_0':
        'quay.io/biocontainers/snp-dists:1.2.0--h577a1d6_0' }"

    input:
    path(alignment)

    output:
    path "*.dists.tsv" , emit: distances
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'cladebreaker'
    """
    snp-dists -q ${args} ${alignment} > ${prefix}.dists.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snp-dists: \$(snp-dists -v 2>&1 | sed 's/^snp-dists //')
    END_VERSIONS
    """
}
