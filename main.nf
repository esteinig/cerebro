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


// include { CultureIdentification } from './lib/production/culture';
// include { BacterialEnrichment } from './lib/production/bacterial';

include { QualityControl; QualityControlNanopore } from './lib/production/quality';
include { PathogenDetection } from './lib/production/pathogen';
include { PanviralEnrichment } from './lib/production/panviral';


include { getReads } from './lib/production/utils'; 
include { getQualityControlDatabases } from './lib/production/utils'; 
include { getPanviralEnrichmentDatabases } from './lib/production/utils'; 
include { getPathogenDetectionDatabases } from './lib/production/utils'; 


process StageInputFiles {

    input:
    path(stageJson)

    output:
    tuple env(PIPELINE), env(SAMPLE_ID), path("*.gz")

    script:

    """
    SAMPLE_ID=\$(cerebro-fs stage --json $stageJson --outdir . --pipeline pipeline.txt)
    PIPELINE=\$(cat pipeline.txt)
    """
}


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

    pipelines = StageInputFiles.out.branch(pipelineSelection);
    
    /* Panviral enrichment */

    def panviralDB = getPanviralEnrichmentDatabases();

    PanviralEnrichment(
        pairedReadsFromStage(pipelines.panviral),
        panviralDB.host, 
        panviralDB.virus, 
        panviralDB.control,
    )

    /* Pathogen detection */
    
    def pathogenDB = getPathogenDetectionDatabases();

    PathogenDetection(
        pairedReadsFromStage(pipelines.pathogen),
        pathogenDB.qualityControl,
    )   

}


workflow pathogen {

    /* Metagenomic diagnostics workflow for pathogen detection and profiling */

    def pathogenDB = getPathogenDetectionDatabases();

    PathogenDetection(
        getReads(
            params.fastqPaired, 
            params.fastqNanopore, 
            params.sampleSheet, 
            params.sampleSheetProduction
        ),
        pathogenDB.qualityControl,
        pathogenDB.taxonomicProfile,
        pathogenDB.metagenomeAssembly,
    )   

}


workflow panviral {

    /* Panviral enrichment probe hybridisation capture panels (e.g. Agilent or Twist) 
       with specific consideration for automated consensus genome assembly and species
       to lineage subtyping schemes.
    */

    def panviralDB = getPanviralEnrichmentDatabases();

    PanviralEnrichment(
        getReads(
            params.fastqPaired, 
            params.fastqNanopore, 
            params.sampleSheet, 
            params.sampleSheetProduction
        ),
        panviralDB.virus, 
        panviralDB.qualityControl,
    )

}


workflow quality {

    /* Read quality control and background coverage + depletion  (host, controls, other) */

    def qualityDB = getQualityControlDatabases();

    if (params.nanopore) {
        QualityControlNanopore(
            getReads(
                null, 
                params.fastqNanopore, 
                params.sampleSheet, 
                params.sampleSheetProduction
            ),
            qualityDB,
        )
    } else {
        QualityControl(
            getReads(
                params.fastqPaired, 
                null, 
                params.sampleSheet, 
                params.sampleSheetProduction
            ),
            qualityDB,
        )
    }
   

}


// workflow bacterial {

//     /* Bacterial enrichment probe hybridisation capture panels (e.g. Agilent or Twist) 
//        with specific consideration for mapped read retention instead of host depletion 
//        and targeted species assembly for antimicrobial resistance or plasmid typing. 
//     */

//     def bacterialDB = getBacterialEnrichmentDatabases();

//     BacterialEnrichment(
//         getReads(
//             params.fastqPaired, 
//             params.fastqNanopore, 
//             params.sampleSheet, 
//             params.sampleSheetProduction
//         ),
//         bacterialDB.host, 
//         bacterialDB.panel, 
//         bacterialDB.control
//     )

// }



// workflow culture {

//     /* Bacterial culture identification (hybrid-) assembly and taxonomic profiling */

//     def pathogenDB = getPathogenDetectionDatabases();

//     CultureIdentification(
//         getReads(
//             params.fastqPaired, 
//             params.fastqNanopore, 
//             params.sampleSheet, 
//             params.sampleSheetProduction
//         ),
//         pathogenDB.qualityControl,
//     )   

// }

