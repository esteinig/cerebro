/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/


include { QualityControl } from "./quality";
include { PipelineConfig } from "./utils";


workflow CultureIdentification {

    take:
        reads
        cultureDatabases
        qualityControlDatabases
        productionConfig
        stagedFileData
    main:

        cerebroWorkflow = "culture-identification"
        workflowStarted = java.time.LocalDateTime.now()

        QualityControl(
            reads,
            qualityControlDatabases
        )


        GenomeAssembly(
            QualityControl.out.reads
        )

        SpeciesTyping(
            GenomeAssembly.out.assembly, 
            cultureDatabases.gtdbDatabase
        )


        VirusRecovery.out.results | map { d -> d[1] } | collect | PanviralTable

        QualityControl.out.results.mix(VirusRecovery.out.results) | groupTuple | ProcessOutput

        PipelineConfig(
            ProcessOutput.out.samples | collect,
            cerebroWorkflow, 
            workflowStarted
        )

}