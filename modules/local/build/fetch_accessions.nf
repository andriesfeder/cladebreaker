process FETCH_ACCESSIONS {
    tag "taxid:${taxid}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::biopython=1.79" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79':
        'quay.io/biocontainers/biopython:1.79' }"

    input:
    val taxid
    val genome_count

    output:
    path "accessions.txt", emit: accessions
    path "versions.yml",   emit: versions

    script:
    def max_records = genome_count > 0 ? genome_count : 100000
    """
    #!/usr/bin/env python3
    from Bio import Entrez
    import sys

    Entrez.email = "cladebreaker@pipeline.org"

    term = (
        "txid${taxid}[Organism] AND latest[filter] "
        "AND all[filter] NOT anomalous[filter]"
    )

    # Search and keep server-side history for large result sets
    handle = Entrez.esearch(db="assembly", term=term, retmax=${max_records}, usehistory="y")
    record = Entrez.read(handle)
    handle.close()

    total    = int(record["Count"])
    webenv   = record["WebEnv"]
    qkey     = record["QueryKey"]
    retmax   = min(total, ${max_records})

    if total == 0:
        sys.exit("ERROR: No assemblies found in NCBI for TaxID ${taxid}")

    print(f"Found {total} assemblies for TaxID ${taxid}; fetching {retmax} accessions ...")

    # Fetch document summaries in batches of 500
    batch   = 500
    accessions = []
    for start in range(0, retmax, batch):
        handle = Entrez.esummary(
            db="assembly",
            retstart=start,
            retmax=min(batch, retmax - start),
            webenv=webenv,
            query_key=qkey,
        )
        summaries = Entrez.read(handle, validate=False)
        handle.close()
        for uid in summaries["DocumentSummarySet"]["DocumentSummary"]:
            acc = uid.get("AssemblyAccession", "")
            if acc.startswith("GCA_") or acc.startswith("GCF_"):
                accessions.append(acc)

    with open("accessions.txt", "w") as fh:
        fh.write("\\n".join(accessions) + "\\n")

    print(f"Wrote {len(accessions)} accessions to accessions.txt")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
