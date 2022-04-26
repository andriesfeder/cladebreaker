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

process SNIPPY_SNIPPY {
    tag "$meta.id"
    label 'process_low'

    //TODO: Fix this
    conda (params.enable_conda ? "bioconda::snippy=4.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_2':
        'quay.io/biocontainers/snippy:4.6.0--hdfd78af_1' }"

    publishDir "${params.outdir}/${meta.id}/snippy", mode: params.publish_dir_mode, overwrite: params.force

    input:

    tuple val(meta), path(snp_in), path(ref)

    output:

    tuple val(meta), path("*.tab")              , optional:true, emit: tab
    tuple val(meta), path("*.csv")              , optional:true, emit: csv
    tuple val(meta), path("*.html")             , optional:true, emit: html
    tuple val(meta), path("*.vcf")              , optional:true, emit: vcf
    tuple val(meta), path("*.bed")              , optional:true, emit: bed
    tuple val(meta), path("*.gff")              , optional:true, emit: gff
    tuple val(meta), path("*.bam")              , optional:true, emit: bam
    tuple val(meta), path("*.bam.bai")          , optional:true, emit: bam_bai
    tuple val(meta), path("*.log")              , emit: log
    tuple val(meta), path("*.aligned.fa")       , optional:true, emit: aln_fa
    tuple val(meta), path("*.consensus.fa")     , optional:true, emit: consensus_fa
    tuple val(meta), path("*.consensus.subs.fa"), optional:true, emit: consensus_subs
    tuple val(meta), path("*.raw.vcf")          , optional:true, emit: raw_vcf
    tuple val(meta), path("*.filt.vcf")         , optional:true, emit: filt_vcf
    tuple val(meta), path("*.vcf.gz")           , optional:true, emit: vcf_gz
    tuple val(meta), path("*.vcf.gz.csi")       , optional:true, emit: vcf_gz_csi
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    def input = ""
    if ( meta.assembly ) {
        input = "--ctgs ${snp_in}"
    }
    else {
        if( meta.single_end ) {
            input = " --R1 ${snp_in[0]}"
        }
        else {
            input = " --R1 ${snp_in[0]} --R2 ${snp_in[1]}"
        }
    }

    """
    snippy \\
        $args \\
        -- ./ \\
        $input \\
        --ref ${ref} \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: snippy: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy //')
    END_VERSIONS
    """
}
