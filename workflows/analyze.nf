/*
========================================================================================
    CLADEBREAKER ANALYZE WORKFLOW

    Runs statistics on an existing phylogenetic tree, without rebuilding it.

    Each test answers a different question, and the default 'auto' runs them all
    and lets the report decide which one leads:

      monophyly         is each labeled group a clade, and if not, what breaks it?
      rosenberg         if it is, could that monophyly be chance? (Rosenberg 2007,
                        with the k-group generalisation of Zhu, Degnan & Steel 2011)
      gsi               how far towards exclusive ancestry has it sorted?
                        (Cummings, Neel & Shaw 2008)
      slatkin_maddison  if it is not a clade, is the residual mixing still
                        significantly less than random? (Slatkin & Maddison 1989)
      snp_separation    how far apart are the groups in SNP distance?

    Everything is cheap next to building the tree, so all selected tests run and
    the decision only chooses which result leads the report.

    Entry point  : nextflow run main.nf --workflow ANALYZE
    Required     : --tree     Newick tree, or a file of trees for an ensemble
                   --groups   Two-column table mapping tip label to group
    Optional     : --tests            auto (default), or a comma-separated subset
                   --alignment        core alignment, enables snp_separation
                   --distances        a precomputed snp-dists matrix instead
                   --gsi_root         as-is (default) | midpoint | outgroup
                   --gsi_outgroup     tip labels for --gsi_root outgroup
                   --gsi_permutations permutations for the significance tests
                   --no-report        skip rendering the PDF report
                   --outdir           output directory (default: ./results)
========================================================================================
*/

nextflow.enable.dsl = 2

include { GSI              } from '../modules/local/gsi/main'
include { MONOPHYLY        } from '../modules/local/analyze/monophyly'
include { ROSENBERG        } from '../modules/local/analyze/rosenberg'
include { SLATKIN_MADDISON } from '../modules/local/analyze/slatkin_maddison'
include { SNP_DISTS        } from '../modules/local/analyze/snp_dists'
include { SNP_SEPARATION   } from '../modules/local/analyze/snp_separation'
include { ANALYZE_REPORT   } from '../modules/local/analyze/report'
include { ANALYZE_PDF      } from '../modules/local/analyze/pdf'

workflow ANALYZE {

    // Declared inside the workflow: NF 26.x forbids top-level statements.
    def ALL_TESTS = ['gsi', 'monophyly', 'rosenberg', 'slatkin_maddison', 'snp_separation']

    if (!params.tree || !params.groups) {
        exit 1, """
        ERROR: --tree and --groups are both required for cladebreaker analyze.

        Usage:
            nextflow run main.nf --workflow ANALYZE \\
                --tree <tree.nwk> \\
                --groups <groups.tsv> \\
                [--tests auto|gsi,monophyly,rosenberg,slatkin_maddison,snp_separation] \\
                [--alignment <core.aln> | --distances <core.dists.tsv>] \\
                [--gsi_root as-is|midpoint|outgroup] \\
                [--gsi_outgroup <tip[,tip...]>] \\
                [--gsi_permutations <N>] \\
                [--outdir <path>]

        The groups file is two columns, tip label and group name:

            tip_label\tgroup
            SAMPLE_1\tquery
            GCA_000012345.1\treference
        """.stripIndent()
    }

    // gsi, monophyly and Rosenberg are all undefined on an unrooted tree, and a
    // mis-rooted tree still yields plausible-looking numbers, so the rooting is
    // always an explicit choice. Slatkin-Maddison alone is rooting independent.
    def rooting = params.gsi_root ?: 'as-is'
    if (!(rooting in ['as-is', 'midpoint', 'outgroup'])) {
        exit 1, "ERROR: --gsi_root must be one of 'as-is', 'midpoint' or 'outgroup' (got '${rooting}')."
    }
    if (rooting == 'outgroup' && !params.gsi_outgroup) {
        exit 1, "ERROR: --gsi_root outgroup requires --gsi_outgroup <tip[,tip...]>."
    }

    def requested = (params.tests ?: 'auto').toString().toLowerCase().tokenize(',')*.trim().findAll { it }
    def unknown   = requested.findAll { !(it in ALL_TESTS + ['auto']) }
    if (unknown) {
        exit 1, "ERROR: --tests does not recognise: ${unknown.join(', ')}.\n" +
                "Valid values are 'auto' or a comma-separated subset of: ${ALL_TESTS.join(', ')}."
    }
    def tests = ('auto' in requested) ? ALL_TESTS : requested.unique()

    if (params.alignment && params.distances) {
        exit 1, "ERROR: pass either --alignment or --distances, not both."
    }
    def wants_snp = 'snp_separation' in tests
    def has_snp_input = params.alignment || params.distances
    if (wants_snp && !has_snp_input && !('auto' in requested)) {
        exit 1, "ERROR: --tests snp_separation needs --alignment <core.aln> or --distances <matrix.tsv>."
    }
    // Under 'auto' the SNP test is simply skipped when no alignment is supplied.
    def run_snp = wants_snp && has_snp_input

    def ch_tree   = file(params.tree,   checkIfExists: true)
    def ch_groups = file(params.groups, checkIfExists: true)

    log.info """\
        ================================================================
         C L A D E B R E A K E R   A N A L Y Z E
        ================================================================
         Tree         : ${params.tree}
         Groups       : ${params.groups}
         Tests        : ${tests.join(', ')}${run_snp ? '' : (wants_snp ? '  (snp_separation skipped: no alignment)' : '')}
         Rooting      : ${rooting}${rooting == 'outgroup' ? " (${params.gsi_outgroup})" : ''}
         Permutations : ${params.gsi_permutations}
         Output dir   : ${params.outdir}
        ================================================================
        """.stripIndent()

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    if ( 'gsi' in tests ) {
        GSI ( ch_tree, ch_groups )
        ch_reports  = ch_reports.mix(GSI.out.json)
        ch_versions = ch_versions.mix(GSI.out.versions)
    }
    if ( 'monophyly' in tests ) {
        MONOPHYLY ( ch_tree, ch_groups )
        ch_reports  = ch_reports.mix(MONOPHYLY.out.json)
        ch_versions = ch_versions.mix(MONOPHYLY.out.versions)
    }
    if ( 'rosenberg' in tests ) {
        ROSENBERG ( ch_tree, ch_groups )
        ch_reports  = ch_reports.mix(ROSENBERG.out.json)
        ch_versions = ch_versions.mix(ROSENBERG.out.versions)
    }
    if ( 'slatkin_maddison' in tests ) {
        SLATKIN_MADDISON ( ch_tree, ch_groups )
        ch_reports  = ch_reports.mix(SLATKIN_MADDISON.out.json)
        ch_versions = ch_versions.mix(SLATKIN_MADDISON.out.versions)
    }
    if ( run_snp ) {
        if ( params.distances ) {
            ch_distances = Channel.fromPath(params.distances, checkIfExists: true)
        } else {
            SNP_DISTS ( file(params.alignment, checkIfExists: true) )
            ch_distances = SNP_DISTS.out.distances
            ch_versions  = ch_versions.mix(SNP_DISTS.out.versions)
        }
        SNP_SEPARATION ( ch_distances, ch_groups )
        ch_reports  = ch_reports.mix(SNP_SEPARATION.out.json)
        ch_versions = ch_versions.mix(SNP_SEPARATION.out.versions)
    }

    ANALYZE_REPORT ( ch_reports.collect() )

    // --no-report and --no_report both work. Nextflow camel-cases a kebab flag,
    // so --no-report arrives as params.noReport, not params.'no-report'.
    def skip_report = params.no_report || (params.containsKey('noReport') && params.noReport)

    // The figure draws the rooted tree MONOPHYLY wrote, which is the tree the
    // statistics actually ran on, so the PDF needs that test to have run.
    if ( !skip_report && 'monophyly' in tests ) {
        ANALYZE_PDF (
            ch_reports.mix(ANALYZE_REPORT.out.json).collect(),
            MONOPHYLY.out.rooted_tree,
            ch_groups
        )
        ch_versions = ch_versions.mix(ANALYZE_PDF.out.versions)
        ANALYZE_PDF.out.pdf
            .map { pdf -> log.info "PDF report written to: ${params.outdir}/analyze/${pdf.name}" }
    }
    else if ( !skip_report ) {
        log.warn "Skipping the PDF report: it draws the rooted tree that the 'monophyly' test writes, which --tests did not select."
    }

    ANALYZE_REPORT.out.report
        .map { report -> log.info "Analyze report written to: ${params.outdir}/analyze/${report.name}" }
}
