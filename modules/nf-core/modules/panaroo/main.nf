process PANAROO {
    label 'process_medium'

    conda (params.enable_conda ? "${moduleDir}/environment.yml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/panaroo:1.6.0--pyhdfd78af_0' :
        'quay.io/biocontainers/panaroo:1.6.0--pyhdfd78af_0' }"

    publishDir { "${params.outdir}/panaroo" }, mode: params.publish_dir_mode, overwrite: params.force

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("results/*")                                       , emit: results
    tuple val(meta), path("results/core_gene_alignment.aln"), optional: true , emit: aln
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Default args include --alignment core so the core gene alignment is produced
    // for downstream phylogenetic inference; override via process.ext.args in modules.config
    def args = task.ext.args ?: '--clean-mode strict --alignment core'
    """
    panaroo \\
        $args \\
        -t $task.cpus \\
        -o results \\
        -i $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panaroo: \$(echo \$(panaroo --version 2>&1) | sed 's/^.*panaroo //')
    END_VERSIONS
    """
}
