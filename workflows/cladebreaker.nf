/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCladebreaker.initialise(params, log)

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
include { INPUT_CHECK } from '../subworkflows/local/input_check'

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

include { WHATSGNU_MAIN                    } from '../modules/local/whatsgnu/main'
include { WHATSGNU_GETGENOMES              } from '../modules/local/whatsgnu/getgenomes'
include { QC_READS                    } from '../modules/local/cladebreaker/qc_reads'

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

    // GATHER_SAMPLES (
    //    INPUT_CHECK.out.reads
        // , params.genome_size
    // )

    //
    // MODULE: qc-reads
    //
    // QC_READS (
    //    GATHER_SAMPLES.out.raw_fastq
    // )

    //
    // MODULE: Run FastQC
    //

    FASTQC (
        INPUT_CHECK.out.reads
    )

    //
    // MODULE: Run Shovill
    //

    SHOVILL (
        INPUT_CHECK.out.reads
    )

    ASSEMBLYSCAN (
        SHOVILL.out.contigs
    )

    //
    //MODULE: Run Prokka
    //
    
    PROKKA (
        SHOVILL.out.contigs.combine(Channel.fromPath( params.proteins )).combine(Channel.fromPath( params.prodigal_tf ))
    )

    //
    //MODULE: Run WhatsGNU
    //

    WHATSGNU_MAIN (
        PROKKA.out.faa.combine(Channel.fromPath( params.db ))
    )

    //
    //MODULE: Run WhatsGNU Gather Genomes
    //

    WHATSGNU_GETGENOMES (
        WHATSGNU_MAIN.out.topgenomes
    )

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    // ch_versions = ch_versions.mix(GATHER_SAMPLES.out.versions.first())
    // ch_versions = ch_versions.mix(QC_READS.out.versions.first())
    ch_versions = ch_versions.mix(SHOVILL.out.versions.first())
    ch_versions = ch_versions.mix(ASSEMBLYSCAN.out.versions.first())
    ch_versions = ch_versions.mix(PROKKA.out.versions.first())
    ch_versions = ch_versions.mix(WHATSGNU_MAIN.out.versions.first())
    ch_versions = ch_versions.mix(WHATSGNU_GETGENOMES.out.versions.first())


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
