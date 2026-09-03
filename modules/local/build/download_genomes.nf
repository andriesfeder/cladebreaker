process DOWNLOAD_GENOMES {
    tag "taxid:${taxid}|${accessions.name}"
    label 'process_low'

    // Cap how many download jobs hit NCBI at once. NCBI throttles aggressive
    // clients with HTTP 503 ("Service unavailable") pages, so this stays low
    // regardless of the executor queueSize. Retry the whole chunk with backoff
    // on a throttle-driven failure (exit 1 is NOT retried by the base config).
    maxForks      params.download_forks ?: 4
    errorStrategy 'retry'
    maxRetries    params.download_max_retries ?: 5

    conda (params.enable_conda ? "bioconda::ncbi-genome-download=0.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.3--pyh864c0ab_1':
        'quay.io/biocontainers/ncbi-genome-download:0.3.3--pyh864c0ab_1' }"

    input:
    path  accessions
    val   taxid

    output:
    path "genomes/*.fna", emit: fna
    path "versions.yml",  emit: versions

    script:
    """
    set -eu
    mkdir -p genomes

    # Back off on retries so we don't immediately re-hit a throttling NCBI.
    if [ ${task.attempt} -gt 1 ]; then
        backoff=\$(( (${task.attempt} - 1) * 120 ))
        echo "Retry attempt ${task.attempt}: sleeping \${backoff}s for NCBI rate limits to recover..."
        sleep \${backoff}
    fi

    total=\$(grep -c . ${accessions} || true)

    # Accessions come in two flavours that live in different NCBI sections:
    #   GCA_* -> genbank      GCF_* -> refseq
    # -s/--section only accepts ONE section per call, so split by prefix and
    # download each from the correct section into the same output folder.
    grep '^GCA_' ${accessions} > gca.txt || true
    grep '^GCF_' ${accessions} > gcf.txt || true

    # ncbi-genome-download verifies md5 and retries internally (-r). Don't abort
    # on a partial failure here; we validate and decide below.
    set +e
    ngd_rc=0
    if [ -s gca.txt ]; then
        ncbi-genome-download \\
            -s genbank \\
            -F fasta \\
            -A gca.txt \\
            --output-folder genomes_raw \\
            --flat-output \\
            -p 1 \\
            -r 10 \\
            bacteria
        ngd_rc=\$(( ngd_rc + \$? ))
    fi
    if [ -s gcf.txt ]; then
        ncbi-genome-download \\
            -s refseq \\
            -F fasta \\
            -A gcf.txt \\
            --output-folder genomes_raw \\
            --flat-output \\
            -p 1 \\
            -r 10 \\
            bacteria
        ngd_rc=\$(( ngd_rc + \$? ))
    fi
    set -e

    # Drop corrupt downloads (e.g. a 503 'Service unavailable' HTML page saved
    # with a .fna.gz name) so they can never poison downstream annotation.
    if [ -d genomes_raw ]; then
        find genomes_raw -name "*.fna.gz" | while read f; do
            gzip -t "\$f" 2>/dev/null || rm -f "\$f"
        done

        # Decompress validated FASTA into a flat genomes/ directory.
        find genomes_raw -name "*.fna.gz" | while read f; do
            gunzip -c "\$f" > "genomes/\$(basename \${f%.gz})"
        done
    fi

    good=\$(ls genomes/*.fna 2>/dev/null | wc -l)
    echo "Downloaded \${good} of \${total} genome FASTA files (ncbi-genome-download rc=\${ngd_rc})."

    # Require most of the chunk to succeed. A large shortfall means NCBI throttled
    # us, so fail to trigger a retry with backoff. A handful of permanently
    # unavailable accessions is tolerated so one dead record can't block the build.
    threshold=\$(( total * 80 / 100 ))
    if [ "\${good}" -lt "\${threshold}" ]; then
        echo "ERROR: only \${good}/\${total} genomes downloaded (< 80%); likely NCBI throttling. Failing for retry." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$( ncbi-genome-download --version )
    END_VERSIONS
    """
}
