/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCladebreaker.initialise(params, log)

errorStrategy 'retry'

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta , params.outdir , params.proteins , params.prodigal_tf , params.db ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { GATHER_GENOMES } from '../subworkflows/local/gather_genomes'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { SHOVILL                     } from '../modules/nf-core/modules/shovill/main'
include { ASSEMBLYSCAN                } from '../modules/nf-core/modules/assemblyscan/main'
include { PROKKA                      } from '../modules/nf-core/modules/prokka/main'
include { ROARY                       } from '../modules/nf-core/modules/roary/main'
include { PIRATE                      } from '../modules/nf-core/modules/pirate/main'
include { RAXMLNG                     } from '../modules/nf-core/modules/raxmlng/main'
include { NCBIGENOMEDOWNLOAD          } from '../modules/nf-core/modules/ncbigenomedownload/main'

include { WHATSGNU_MAIN               } from '../modules/local/whatsgnu/main'
include { WHATSGNU_GETGENOMES         } from '../modules/local/whatsgnu/getgenomes'
include { QC_READS                    } from '../modules/local/cladebreaker/qc_reads'
include { SNIPPY                      } from '../modules/local/snippy/snippy'
include {SNIPPY_CORE                  } from '../modules/local/snippy/snippycore'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CLADEBREAKER {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // MODULE: Run FastQC
    // Need to figure out the genome_size, setup_datasets thing to eventually
    // replace params.genome_size

    //
    // MODULE: Run FastQC
    //
    assemblyscan_input = Channel.empty()
    FASTQC (
        INPUT_CHECK.out.reads
    )

    //
    // MODULE: Run Shovill
    //

    SHOVILL (
        INPUT_CHECK.out.reads
    )
    assemblyscan_input = assemblyscan_input.mix(SHOVILL.out.contigs)

    assemblyscan_input = assemblyscan_input.mix(INPUT_CHECK.out.assemblies)
    ASSEMBLYSCAN (
        assemblyscan_input
    )

    //
    //MODULE: Run Prokka
    //
    prokka_input = Channel.empty()
    prokka_input = prokka_input.mix(INPUT_CHECK.out.assemblies)
    prokka_input = prokka_input.mix(SHOVILL.out.contigs)
    // prokka_input = prokka_input.combine(Channel.fromPath( params.proteins )).combine(Channel.fromPath( params.prodigal_tf ))

    PROKKA (
        prokka_input.combine(Channel.fromPath( params.proteins )).combine(Channel.fromPath( params.prodigal_tf ))
    )

    //
    //MODULE: Run WhatsGNU
    //

    WHATSGNU_MAIN (
        PROKKA.out.faa.combine(Channel.fromPath( params.db ))
    )

    //
    //SUBWORKFLOW: Run NCBI Genome Download and Prokka for each genbank genome
    //

    GATHER_GENOMES (
        WHATSGNU_MAIN.out.gca_list
    )

    //
    //MODULE: Run Roary
    //
    if ( params.ref == null) {
        roary_input = Channel.empty()
        roary_input = roary_input.mix(GATHER_GENOMES.out.prokka_gff.last())
        roary_input = roary_input.combine(Channel.fromPath("${workflow.workDir}/tmp/gff/"))

        ROARY (
            roary_input
        )
        if ( params.run_raxml ) {
            RAXMLNG (
                ROARY.out.aln
            )
        }
    }
    else {
        snippy_input = Channel.empty()
        snippy_input = snippy_input.mix(INPUT_CHECK.out.reads).mix(INPUT_CHECK.out.assemblies)
        snippy_input = snippy_input.mix(GATHER_GENOMES.out.ncbi)
        snippy_input = snippy_input.combine(Channel.fromPath( params.ref ))
        SNIPPY (
            snippy_input
        )
        snippy_core = Channel.empty()
        snippy_core = SNIPPY.out.snippy_dir.collect()
        SNIPPY_CORE (
            snippy_core
        )
        if ( params.run_raxml ) {
            RAXMLNG (
                SNIPPY_CORE.out.full_aln
            )
        }
    }

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_versions = ch_versions.mix(SHOVILL.out.versions.first())
    ch_versions = ch_versions.mix(ASSEMBLYSCAN.out.versions.first())
    ch_versions = ch_versions.mix(PROKKA.out.versions.first())
    ch_versions = ch_versions.mix(WHATSGNU_MAIN.out.versions.first())
    ch_versions = ch_versions.mix(GATHER_GENOMES.out.versions.first())
    if ( params.ref == null) {
        ch_versions = ch_versions.mix(ROARY.out.versions)
    }
    else {
        ch_versions = ch_versions.mix(SNIPPY.out.versions)
        ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions)
    }
    if ( params.run_raxml ){
        ch_versions = ch_versions.mix(RAXMLNG.out.versions)
    }
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCladebreaker.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
