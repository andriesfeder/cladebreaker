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

process SNIPPY {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::snippy=4.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_2':
        'quay.io/biocontainers/snippy:4.6.0--hdfd78af_1' }"

    publishDir "${params.outdir}/${meta.id}/snippy/", mode: params.publish_dir_mode, overwrite: params.force

    input:

    tuple val(meta), path(snp_in), path(ref)

    output:

    tuple val(meta), path("${meta.id}/*.tab")              , optional:true, emit: tab
    tuple val(meta), path("${meta.id}/*.csv")              , optional:true, emit: csv
    tuple val(meta), path("${meta.id}/*.html")             , optional:true, emit: html
    tuple val(meta), path("${meta.id}/*.vcf")              , optional:true, emit: vcf
    tuple val(meta), path("${meta.id}/*.bed")              , optional:true, emit: bed
    tuple val(meta), path("${meta.id}/*.gff")              , optional:true, emit: gff
    tuple val(meta), path("${meta.id}/*.bam")              , optional:true, emit: bam
    tuple val(meta), path("${meta.id}/*.bam.bai")          , optional:true, emit: bam_bai
    tuple val(meta), path("${meta.id}/*.log")              , emit: log
    tuple val(meta), path("${meta.id}/*.aligned.fa")       , optional:true, emit: aln_fa
    tuple val(meta), path("${meta.id}/*.consensus.fa")     , optional:true, emit: consensus_fa
    tuple val(meta), path("${meta.id}/*.consensus.subs.fa"), optional:true, emit: consensus_subs
    tuple val(meta), path("${meta.id}/*.raw.vcf")          , optional:true, emit: raw_vcf
    tuple val(meta), path("${meta.id}/*.filt.vcf")         , optional:true, emit: filt_vcf
    tuple val(meta), path("${meta.id}/*.vcf.gz")           , optional:true, emit: vcf_gz
    tuple val(meta), path("${meta.id}/*.vcf.gz.csi")       , optional:true, emit: vcf_gz_csi
    path "${meta.id}/"                                     , emit: snippy_dir
    path "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
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
        --outdir ${meta.id}/ \\
        $input \\
        --reference ${ref} \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: \$(echo \$(snippy --version 2>&1) | sed 's/^.*snippy //')
    END_VERSIONS
    """
}
