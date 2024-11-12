/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/


include { QualityControl } from "./quality";
include { VirusRecovery; ProcessOutput } from "../processes/panviral";


workflow PanviralEnrichment {

    take:
        reads
        virusDatabase
        qualityControlDatabases
    main:

        QualityControl(
            reads,
            qualityControlDatabases
        )

        VirusRecovery(
            QualityControl.out.reads, 
            virusDatabase, 
            params.panviralEnrichment.virusAligner,
            params.panviralEnrichment.vircovArgs
        )

        results = QualityControl.out.results.mix(
            VirusRecovery.out.results
        )

        results | groupTuple | ProcessOutput

}