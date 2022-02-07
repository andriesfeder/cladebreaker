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
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { CLADEBREAKER } from './workflows/cladebreaker'

//
// WORKFLOW: Run main nf-core/cladebreaker analysis pipeline
//
workflow NFCORE_CLADEBREAKER {
    CLADEBREAKER ()
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
    NFCORE_CLADEBREAKER ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
