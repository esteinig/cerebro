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


params.outdir = "workflow"

params.pipeline = [
    id: "6114563e-1052-4d5b-930a-1de496f2b65e",
    name: "TWIST-v1",
    location: "DGX",
    stage: "6114563e-1052-4d5b-930a-1de496f2b65e",
]

params.team = "CNS"
params.stageDirectory = "test"
params.deleteStaged = true
params.stageInterval = 10

include { StageInputFiles } from './lib/production/utils';
include { PanviralEnrichment } from './lib/production/panviral';

include { stageCerebroSamples } from './lib/production/utils';
include { getPanviralEnrichmentHostDatabase } from './lib/production/utils'; 
include { getPanviralEnrichmentVirusDatabase } from './lib/production/utils'; 



def pipelineSelection = branchCriteria {
    panviral: it[0] == 'panviral-enrichment'
    pathogen: it[0] == 'pathogen-detection'
    culture:  it[0] == 'culture-identification'
}

def pairedReadsFromStage(staged) {
    return tuple(
        staged[1], 
        staged[2], 
        staged[3], 
        staged[4][0], 
        staged[4][1]
    )
}

workflow production {

    staged = Channel.watchPath("$params.stageDirectory/*.json") | StageInputFiles 

    stageCerebroSamples(
        params.team, 
        params.pipeline.id, 
        params.stageDirectory, 
        params.deleteStaged, 
        params.stageInterval
    );

    /* Pipeline branches and launch */

    pipelines = staged.branch(pipelineSelection)
        
    PanviralEnrichment(
        pipelines.panviral.map { pairedReadsFromStage(it) },
        getPanviralEnrichmentHostDatabase(),
        getPanviralEnrichmentVirusDatabase()
    )
    
}