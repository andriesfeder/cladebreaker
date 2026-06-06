process BAKTA {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bakta" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.9.4--pyhdfd78af_0':
        'quay.io/biocontainers/bakta:1.9.4--pyhdfd78af_0' }"

    publishDir { "${params.outdir}/${meta.id}/" }, mode: params.publish_dir_mode, overwrite: params.force

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("annotation/${meta.id}.faa"), emit: faa
    tuple val(meta), path("annotation/${meta.id}.gff3"), emit: gff
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    bakta \\
        --db ${db} \\
        --output annotation \\
        --prefix ${meta.id} \\
        --threads ${task.cpus} \\
        ${args} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version 2>&1) | sed 's/bakta //')
    END_VERSIONS
    """
}
