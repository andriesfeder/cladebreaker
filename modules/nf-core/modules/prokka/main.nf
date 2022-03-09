process PROKKA {
    tag "$meta.id"
    label 'process_low'

    //TODO: Figure out conda build issue.
    // conda (params.enable_conda ? "bioconda::prokka=1.14.6" : null)
    // conda (params.enable_conda ? "${baseDir}/modules/nf-core/modules/prokka/prokka-env.yaml" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0' :
        'quay.io/biocontainers/prokka:1.14.6--pl526_0' }"

    publishDir "${params.outdir}/${meta.id}/annotation", mode: params.publish_dir_mode, overwrite: params.force

    input:
    tuple val(meta), path(fasta), path(proteins), path(prodigal_tf)

    output:
    tuple val(meta), path("*.gff"), emit: gff
    tuple val(meta), path("*.gbk"), emit: gbk
    tuple val(meta), path("*.fna"), emit: fna
    tuple val(meta), path("*.faa"), emit: faa
    tuple val(meta), path("*.ffn"), emit: ffn
    tuple val(meta), path("*.sqn"), emit: sqn
    tuple val(meta), path("*.fsa"), emit: fsa
    tuple val(meta), path("*.tbl"), emit: tbl
    tuple val(meta), path("*.err"), emit: err
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    prokka_prodigal = ""
    if (prodigal_tf.getName() != 'EMPTY_TF' && !params.skip_prodigal_tf) {
        prokka_prodigal = "--prodigaltf ${prodigal_tf}"
    }

    prokka_proteins = ""
    if (proteins.getName() != 'EMPTY_PROTEINS') {
        prokka_proteins = "--proteins ${proteins}"
        proteins_name = prokka_proteins.getName()
        if(proteins_name.contains("-")) {
            genus = proteins_name.split('-')[0].capitalize()
            species = proteins_name.split('-')[1]
        } else {
            genus = proteins_name.capitalize()
            species = "spp."
        }
    }

    // def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    // def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    //${params.outdir}/${meta.id}/$fasta
    """
    prokka \\
        $args \\
        --cpus $task.cpus \\
        --prefix $prefix \\
        --outdir ./ \\
        --force \\
        $prokka_proteins \\
        $prokka_prodigal \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
