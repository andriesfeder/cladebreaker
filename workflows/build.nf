/*
========================================================================================
    CLADEBREAKER BUILD WORKFLOW

    Downloads reference genomes from NCBI for a given TaxID, annotates them
    with Prokka, and builds a WhatsGNU database using WhatsGNU_database_customizer
    and WhatsGNU_main.

    Entry point  : nextflow run main.nf -entry CLADEBREAKER_BUILD
    Required     : --taxid       NCBI TaxID of the target species
    Optional     : --genome_count  Max genomes to use (default: all available)
                   --db_mode       'basic' or 'ortholog'  (default: basic)
                   --outdir        Output directory (default: ./cladebreaker_databases)
========================================================================================
*/

nextflow.enable.dsl = 2

include { FETCH_ACCESSIONS  } from '../modules/local/build/fetch_accessions'
include { DOWNLOAD_GENOMES  } from '../modules/local/build/download_genomes'
include { PROKKA            } from '../modules/nf-core/modules/prokka/main'
include { WHATSGNU_BUILD_DB } from '../modules/local/whatsgnu/database_customizer'

workflow CLADEBREAKER_BUILD {

    if (!params.taxid) {
        exit 1, """
        ERROR: --taxid is required for cladebreaker build.

        Usage:
            nextflow run main.nf -entry CLADEBREAKER_BUILD \\
                --taxid <NCBI_TAXID> \\
                [--genome_count <N>] \\
                [--db_mode basic|ortholog] \\
                [--outdir <path>]

        Tip: Use cladebreaker_build.py for an interactive guided experience.
        """.stripIndent()
    }

    def genome_count = params.genome_count ?: 0   // 0 signals "all available"
    def db_mode      = params.db_mode      ?: "basic"

    log.info """\
        ================================================================
         C L A D E B R E A K E R   B U I L D
        ================================================================
         TaxID        : ${params.taxid}
         Genome count : ${genome_count > 0 ? genome_count : 'all available'}
         DB mode      : ${db_mode}
         Output dir   : ${params.outdir}
        ================================================================
        """.stripIndent()

    //
    // Fetch genome accession list from NCBI Assembly for this TaxID
    //
    FETCH_ACCESSIONS(
        params.taxid,
        genome_count
    )

    //
    // Download all accessions in one ncbi-genome-download call
    //
    DOWNLOAD_GENOMES(
        FETCH_ACCESSIONS.out.accessions,
        params.taxid
    )

    //
    // Fan out: one channel item per downloaded FNA file → annotate with Prokka
    //
    ch_for_prokka = DOWNLOAD_GENOMES.out.fna
        .flatten()
        .map { fna ->
            def meta = [id: fna.baseName.replaceAll(/_genomic$/, '')]
            tuple(
                meta,
                fna,
                file("${workflow.projectDir}/data/EMPTY_PROTEINS"),
                file("${workflow.projectDir}/data/EMPTY_TF")
            )
        }

    PROKKA(ch_for_prokka)

    //
    // Collect all Prokka FAA files → build WhatsGNU database
    //
    ch_faa = PROKKA.out.faa
        .map { meta, faa -> faa }
        .collect()

    WHATSGNU_BUILD_DB(
        ch_faa,
        params.taxid,
        db_mode
    )

    WHATSGNU_BUILD_DB.out.database
        .map { db -> log.info "WhatsGNU database ready: ${db}" }
}
