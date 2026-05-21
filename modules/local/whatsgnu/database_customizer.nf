process WHATSGNU_BUILD_DB {
    tag "taxid:${taxid}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::whatsgnu=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/whatsgnu:1.3--hdfd78af_0':
        'quay.io/biocontainers/whatsgnu:1.3--hdfd78af_0' }"

    publishDir "${params.outdir}/whatsgnu_db", mode: params.publish_dir_mode, overwrite: params.force

    input:
    path   faas        // collected list of Prokka FAA files
    val    taxid
    val    db_mode     // 'basic' or 'ortholog'

    output:
    path "*.pickle",  emit: database
    path "*.txt",     emit: database_txt, optional: true
    path "versions.yml", emit: versions

    script:
    def prefix = "${taxid}_${db_mode}"
    """
    # ----------------------------------------------------------------
    # Collect all FAA files into a single input directory
    # ----------------------------------------------------------------
    mkdir -p faa_input
    for faa in ${faas}; do
        cp \$faa faa_input/
    done

    faa_count=\$(ls faa_input/*.faa 2>/dev/null | wc -l)
    echo "Building WhatsGNU database from \${faa_count} FAA files ..."

    # ----------------------------------------------------------------
    # Step 1: Add genome-name prefix to each protein (customizer)
    # ----------------------------------------------------------------
    WhatsGNU_database_customizer.py \\
        ${prefix} \\
        faa_input/ \\
        -p \\
        -c

    # ----------------------------------------------------------------
    # Step 2: Compress into a WhatsGNU database
    # Locate the concatenated FAA produced by the customizer
    # ----------------------------------------------------------------
    CONCAT=\$(find ${prefix}/ -name "*.faa" | head -1)

    if [ -z "\$CONCAT" ]; then
        echo "ERROR: Could not find concatenated FAA from WhatsGNU_database_customizer.py"
        exit 1
    fi

    echo "Building database from \$CONCAT ..."

    WhatsGNU_main.py \\
        \$CONCAT \\
        -m \$CONCAT \\
        -a \\
        -p ${prefix} \\
        -o WhatsGNU_build_output \\
        --force

    # ----------------------------------------------------------------
    # Locate and rename the output pickle
    # ----------------------------------------------------------------
    PICKLE=\$(find . -maxdepth 2 -name "*.pickle" | head -1)
    if [ -z "\$PICKLE" ]; then
        echo "ERROR: No pickle database file was created."
        exit 1
    fi

    DEST="${prefix}.pickle"
    if [ "\$PICKLE" != "./\$DEST" ]; then
        mv "\$PICKLE" "\$DEST"
    fi

    echo "Database created: \$DEST"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        WhatsGNU: \$(echo \$(WhatsGNU_main.py --version 2>&1) | sed 's/^.*WhatsGNU //')
    END_VERSIONS
    """
}
