process DOWNLOAD_GENOMES {
    tag "taxid:${taxid}"
    label 'process_low'
    maxRetries 3

    conda (params.enable_conda ? "bioconda::ncbi-genome-download=0.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.0--pyh864c0ab_1':
        'quay.io/biocontainers/ncbi-genome-download:0.3.0--pyh864c0ab_1' }"

    input:
    path  accessions
    val   taxid

    output:
    path "genomes/*.fna", emit: fna
    path "versions.yml",  emit: versions

    script:
    """
    mkdir -p genomes

    ncbi-genome-download \\
        -s genbank \\
        -F fasta \\
        -A ${accessions} \\
        --output-folder genomes_raw \\
        --flat-output \\
        -r 3 \\
        bacteria

    # Decompress and move FNA files to a flat genomes/ directory
    find genomes_raw -name "*.fna.gz" | while read f; do
        gunzip -c "\$f" > "genomes/\$(basename \${f%.gz})"
    done

    count=\$(ls genomes/*.fna 2>/dev/null | wc -l)
    echo "Downloaded \${count} genome FASTA files."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$( ncbi-genome-download --version )
    END_VERSIONS
    """
}
