/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
// WorkflowCladebreaker.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta]
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

//include { get_resources; get_schemas; print_efficiency } from '../lib/nf/functions'
include { create_input_channel; setup_datasets } from '../lib/nf/bactopia'
include { get_resources; print_efficiency } from '../lib/nf/functions'
RESOURCES = get_resources(workflow.profile, params.max_memory, params.max_cpus)

//SCHEMAS = get_schemas()
WorkflowCladebreaker.initialise(params, log)
runtype = WorkflowCladebreaker.initialise(params, log)
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
// include { WHATSGNU } from '../subworkflows/local/whatsgnu'
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
include { PROKKA                      } from '../modules/nf-core/modules/prokka/main'
include { ROARY                       } from '../modules/nf-core/modules/roary/main'
include { PIRATE                      } from '../modules/nf-core/modules/pirate/main'
include { RAXMLNG                     } from '../modules/nf-core/modules/raxmlng/main'

include { GATHER_SAMPLES              } from '../modules/local/cladebreaker/gather_samples/main'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CLADEBREAKER {

    datasets = setup_datasets()
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run gather_samples
    //
    GATHER_SAMPLES ( 
        INPUT_CHECK.out.reads
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        // GATHER_SAMPLES.out.raw_fastq
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //CUSTOM_DUMPSOFTWAREVERSIONS (
    //    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    //)

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCladebreaker.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    //
    // MODULE: Shovill
    //
    ch_shovill_files = Channel.empty()

    //SHOVILL (
    //    INPUT_CHECK.out.reads
    //)
    //ch_versions = ch_versions.mix(SHOVILL.out.versions.first())

    //
    // MODULE: Prokka
    //

    //PROKKA (
    //    SHOVILL.out.contigs
    //)
    //ch_versions = ch_versions.mix(PROKKA.out.versions.first())

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
