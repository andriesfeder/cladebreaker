#!/usr/bin/env nextflow
/*
========================================================================================
    cladebreaker
========================================================================================
    Github : https://github.com/andriesfeder/cladebreaker
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
include { BUILD              } from './workflows/build'

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
    CLADEBREAKER ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
