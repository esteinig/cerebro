#!/usr/bin/env nextflow

/* 
vim: syntax=groovy
-*- mode: groovy;-*-

Cerebro: metagenomic and -transcriptomic host-pathogen diagnostics for clinical production environments

=================
Production module
=================

Production is integrated with the Cerebro application stack and command-line interface. It is primarily designed to
operate on a local server in a laboratory or clinical setting, but can be configured as described in the documentation.

Production pipeline operate as follows:





*/


include { StageInputFiles } from './lib/production/utils';
include { PanviralEnrichment } from './lib/production/panviral';

include { getPanviralEnrichmentDatabases } from './lib/production/utils'; 

def pipelineSelection = branchCriteria {
    panviral: it[0] == 'panviral-enrichment'
    pathogen: it[0] == 'pathogen-detection'
    culture:  it[0] == 'culture-identification'
}

def pairedReadsFromStage(channel) {
    return channel.map { tuple(
        it[1], it[2][0], it[2][1]
    ) }
}

workflow production {

    /* Staged sample files provided by tower */

    Channel.fromPath("$params.stageDirectory/*.json") | StageInputFiles

    /* Pipeline selection */
   
    pipelines = StageInputFiles.out.branch(pipelineSelection)
    
    /* Panviral enrichment */

    def panviralDB = getPanviralEnrichmentDatabases()

    PanviralEnrichment(
        pairedReadsFromStage(pipelines.panviral),
        panviralDB.host, 
        panviralDB.virus, 
        panviralDB.control,
    )

    /* Pathogen detection */
    
    

    // def database = getPathogenDetectionDatabases()
    
    // PathogenDetection(

    // )
    
    

}