process WHATSGNU_MAIN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::whatsgnu=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatsgnu:1.3--hdfd78af_0':
        'quay.io/biocontainers/whatsgnu:1.3--hdfd78af_0' }"

    publishDir "${params.outdir}/${meta.id}/WhatsGNU", mode: params.publish_dir_mode, overwrite: params.force

    input:
    tuple val(meta), path(faa), path(database)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("WhatsGNU_Report/*_WhatsGNU_report.txt")  , emit: report
    tuple val(meta), path("WhatsGNU_Report/*_topgenomes.txt")       , emit: topgenomes
    path("WhatsGNU_Report/${meta.id}_gca_list.txt"), emit: gca_list
    tuple val(meta), path("WhatsGNU_Report/*.log")                  , emit: log
    path "versions.yml"                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    mode = ""
    if ( params.o ) {
        mode = "-dm ortholog"
    }
    if ( params.b ) {
        mode = "-dm basic"
    }

    """
    WhatsGNU_main.py \\
        ${args} \\
        -d ${database} \\
        ${mode} \\
        -t \\
        -o WhatsGNU_Report \\
        --topgenomes_count ${params.topgenomes_count} \\
        --force \\
        ${faa}

    list_fixer.py WhatsGNU_Report/*_topgenomes.txt WhatsGNU_Report/${meta.id}_gca_list.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        WhatsGNU: \$(echo \$(WhatsGNU_main.py --version 2>&1) | sed 's/^.*WhatsGNU //')
    END_VERSIONS
    """
}
