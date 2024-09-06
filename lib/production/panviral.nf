/* ==================================
 * PANVIRAL ENRICHMENT FOR PRODUCTION
 * ==================================
*/


include { 
    QualityControl
} from "./quality";


include { 
    VirusRecovery
    ProcessOutput
} from "./panviral";


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
            params.panviralEnrichment.virusAligner
        )

        results = QualityControl.out.results.mix(
            VirusRecovery.out.results
        )

        results | groupTuple | ProcessOutput

}