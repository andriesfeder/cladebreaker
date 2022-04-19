process RAXMLNG {
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::raxml-ng=1.0.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/raxml-ng:1.0.3--h32fcf60_0' :
        'quay.io/biocontainers/raxml-ng:1.0.3--h32fcf60_0' }"

    
    publishDir "${params.outdir}/raxml", mode: params.publish_dir_mode, overwrite: params.force

    input:
    tuple val(meta), path(alignment)


    output:
    path "*.raxml.bestTree", emit: phylogeny
    path "*.raxml.support" , optional:true, emit: phylogeny_bootstrapped
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    raxml-ng \\
        $args \\
        --model GTR+G \\
        --msa $alignment \\
        --threads $task.cpus \\
        --prefix output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raxmlng: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """
}
