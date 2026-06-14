/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/


include { QualityControl } from "./quality";
include { PipelineConfig; ToolVersions; Checksums } from "./utils";
include { VirusRecovery; PanviralTable; ProcessOutput  } from "../processes/panviral";


workflow PanviralEnrichment {

    take:
        reads
        panviralDatabases
        qualityControlDatabases
        productionConfig
        stagedFileData
    main:

        cerebroWorkflow = "panviral-enrichment"
        workflowStarted = java.time.LocalDateTime.now()

        QualityControl(
            reads,
            qualityControlDatabases
        )

        VirusRecovery(
            QualityControl.out.reads, 
            panviralDatabases.virusDatabase, 
            params.panviralEnrichment.virusAligner,
            params.panviralEnrichment.vircovArgs
        )


        VirusRecovery.out.results | map { d -> d[1] } | collect | PanviralTable

        QualityControl.out.results.mix(VirusRecovery.out.results) | groupTuple | ProcessOutput

        PipelineConfig(
            ProcessOutput.out.samples | collect,
            cerebroWorkflow, 
            workflowStarted
        )

        // Lifecycle metadata sidecars (intermission-3; additive, no output change)
        if (params.cerebro.emitToolVersions) { ToolVersions() }
        if (params.cerebro.emitChecksums)    { Checksums(PipelineConfig.out.config) }

}