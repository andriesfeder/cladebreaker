#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/cladebreaker
========================================================================================
    Github : https://github.com/nf-core/cladebreaker
    Website: https://nf-co.re/cladebreaker
    Slack  : https://nfcore.slack.com/channels/cladebreaker
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { CLADEBREAKER       } from './workflows/cladebreaker'
include { CLADEBREAKER_BUILD } from './workflows/build'

//
// WORKFLOW: Run main cladebreaker analysis pipeline
//
workflow NFCORE_CLADEBREAKER {
    CLADEBREAKER ()
}

//
// WORKFLOW: Build a WhatsGNU database for a given NCBI TaxID
//
workflow NFCORE_BUILD {
    CLADEBREAKER_BUILD ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    WorkflowMain.initialise(workflow, params, log)
    NFCORE_CLADEBREAKER ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
