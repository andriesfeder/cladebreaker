// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

process SNIPPY_CORE {
    // tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::snippy=4.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_2':
        'quay.io/biocontainers/snippy:4.6.0--hdfd78af_1' }"

    publishDir "${params.outdir}/snippy-core", mode: params.publish_dir_mode, overwrite: params.force


    input:

    path(paths)

    output:
    // tuple val(meta), path("*.aln")                       , emit: aln
    path("*.aln")                       , emit: aln
    path("*.full.aln")                  , emit: full_aln
    path("*.tab")                       , emit: tab
    path("*.vcf")                       , emit: vcf
    path("*.txt")                       , emit: txt
    path("*.self_mask.bed")             , optional: true, emit: bed
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def paths_in = ""
    for(i in paths){
        paths_in = paths_in + " " + i
    }

    //TODO: deal with meta for snippycore
    def meta = [:]
    meta.id     = "all_samples"
    meta.single_end   = null
    meta.assembly     = null

    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    """
    snippy-core \\
        $args \\
        --ref "${params.ref}"\\
        $paths_in

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy-core: \$(echo \$(snippy-core --version 2>&1) | sed 's/^.*snippy-core //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
