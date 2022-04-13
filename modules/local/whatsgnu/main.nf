// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process WHATSGNU_MAIN {
    tag "$meta.id"
    label 'process_low'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
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

    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO: change --topgenomes_count to ${params.topgenomes_count}
    """
    WhatsGNU_main.py \\
        ${args} \\
        -d ${database} \\
        ${mode} \\
        -t \\
        -o WhatsGNU_Report \\
        --topgenomes_count 3 \\
        --force \\
        ${faa}

    list_fixer.py WhatsGNU_Report/*_topgenomes.txt WhatsGNU_Report/${meta.id}_gca_list.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        WhatsGNU: \$(echo \$(WhatsGNU_main.py --version 2>&1) | sed 's/^.*WhatsGNU //')
    END_VERSIONS
    """
}
