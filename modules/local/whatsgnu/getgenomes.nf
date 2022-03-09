// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string. 
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.

process WHATSGNU_GETGENOMES {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::whatsgnu=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatsgnu:1.3--hdfd78af_0':
        'quay.io/biocontainers/whatsgnu:1.3--hdfd78af_0' }"

    publishDir "${params.outdir}/${meta.id}/WhatsGNU", mode: params.publish_dir_mode, overwrite: params.force


    input:
    
    tuple val(meta), path(topgenomes)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("genbank_genomes/*.fna")  , emit: genbank_genomes
    tuple val(meta), path("${meta.id}_gca_list.txt"), emit: cleanList
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    
    
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    
    """
    WhatsGNU_get_GenBank_genomes.py \\
        -c ${meta.id}_gca_list.txt \\
        genbank_genomes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        WhatsGNU: \$(echo \$(WhatsGNU_main.py --version 2>&1) | sed 's/^.*WhatsGNU //')
    END_VERSIONS 
    """
}
