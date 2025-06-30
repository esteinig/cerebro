/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/


include { QualityControl } from "./quality";
include { PipelineConfig } from "./utils";
include { VirusRecovery; PanviralTable; ProcessOutput; UploadOutput  } from "../processes/panviral";


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


        if (params.cerebroProduction.enabled && params.cerebroProduction.fsConfig.enabled) {
            
            // Collect and upload all output files to CerebroFS

        }

        if (params.cerebroProduction.enabled && params.cerebroProduction.uploadConfig.enabled) {
            
            UploadOutput(
                ProcessOutput.out.results.mix(stagedFileData) | groupTuple | map { d -> d.flatten() }, 
                panviralDatabases.taxonomy, 
                PipelineConfig.out.config
            )

        }

}