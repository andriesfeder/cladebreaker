/*
========================================================================================
    CLADEBREAKER ANALYZE WORKFLOW

    Runs statistics on an existing phylogenetic tree, without rebuilding it.

    Currently the genealogical sorting index (gsi) of Cummings, Neel & Shaw
    (2008), which quantifies how far each labeled group on a rooted tree has
    sorted towards exclusive ancestry, from 0 (mixed) to 1 (monophyletic).

    Entry point  : nextflow run main.nf --workflow ANALYZE
    Required     : --tree     Newick tree, or a file of trees for an ensemble
                   --groups   Two-column table mapping tip label to group
    Optional     : --gsi_root         as-is (default) | midpoint | outgroup
                   --gsi_outgroup     tip labels for --gsi_root outgroup
                   --gsi_permutations Permutations for the significance test
                   --outdir           Output directory (default: ./results)
========================================================================================
*/

nextflow.enable.dsl = 2

include { GSI } from '../modules/local/gsi/main'

workflow ANALYZE {

    if (!params.tree || !params.groups) {
        exit 1, """
        ERROR: --tree and --groups are both required for cladebreaker analyze.

        Usage:
            nextflow run main.nf --workflow ANALYZE \\
                --tree <tree.nwk> \\
                --groups <groups.tsv> \\
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

    // gsi is undefined on an unrooted tree, and a mis-rooted tree still yields
    // plausible-looking numbers, so the rooting is always an explicit choice.
    def rooting = params.gsi_root ?: 'as-is'
    if (!(rooting in ['as-is', 'midpoint', 'outgroup'])) {
        exit 1, "ERROR: --gsi_root must be one of 'as-is', 'midpoint' or 'outgroup' (got '${rooting}')."
    }
    if (rooting == 'outgroup' && !params.gsi_outgroup) {
        exit 1, "ERROR: --gsi_root outgroup requires --gsi_outgroup <tip[,tip...]>."
    }

    def ch_tree   = file(params.tree,   checkIfExists: true)
    def ch_groups = file(params.groups, checkIfExists: true)

    log.info """\
        ================================================================
         C L A D E B R E A K E R   A N A L Y Z E
        ================================================================
         Tree         : ${params.tree}
         Groups       : ${params.groups}
         Rooting      : ${rooting}${rooting == 'outgroup' ? " (${params.gsi_outgroup})" : ''}
         Permutations : ${params.gsi_permutations}
         Output dir   : ${params.outdir}
        ================================================================
        """.stripIndent()

    GSI (
        ch_tree,
        ch_groups
    )

    GSI.out.results
        .map { results -> log.info "Genealogical sorting index written to: ${params.outdir}/gsi/${results.name}" }
}
