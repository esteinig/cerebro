/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/


include { QualityControl } from "./quality";
include { VirusRecovery; PanviralTable; ProcessOutput; UploadOutput  } from "../processes/panviral";


workflow PanviralEnrichment {

    take:
        reads
        panviralDatabases
        qualityControlDatabases
    main:

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

        json = QualityControl.out.results.mix(VirusRecovery.out.results) | groupTuple | ProcessOutput

        if (params.cerebroProduction.enabled) {

            UploadOutput(
                ProcessOutput.out.results, 
                panviralDatabases.taxonomy, 
                "CNS", 
                "CNS", 
                "Default"
            )

        }

}