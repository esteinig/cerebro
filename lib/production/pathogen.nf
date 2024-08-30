/* ==================================
 * PATHOGEN DETECTION FOR PRODUCTION
 * ==================================
*/

include { 
    ReadScan
    SyntheticControls
    ReadQuality
    Deduplication
    HostDepletion
    InternalControls
    BackgroundDepletion
} from "../processes/pathogen"


workflow PathogenDetection {

    take:
        reads
        qualityControlDatabases
    main:

        QualityControl(reads, qualityControlDatabases) | view

}

workflow QualityControl {
    take:
        reads
        databases
    main:

        qualityControlParams = params.pathogenDetection.qualityControl

        // total read counts and table before depletions

        reads = ReadScan(reads).out.reads

        if ( // first for synthetic controls for biomass estimates
            params.pathogenDetection
                .qualityControl
                .syntheticControls
        ) {
            reads = SyntheticControls(
                reads,
                databases.syntheticControls,
                qualityControlParams.syntheticControlsAligner
            ).out.reads
        }

        if ( // second for read quality control and trimming
            params.pathogenDetection
                .qualityControl
                .readQuality
        ) {
            reads = ReadQuality(reads).out.reads
        }

        if (  // third for read head + umi deduplication
            params.pathogenDetection
                .qualityControl
                .readDeduplication
        ) {
            reads = Deduplication(
                reads,
                qualityControlParams.readDeduplicationUmi,
                qualityControlParams.readDeduplicationHead,
                qualityControlParams.readDeduplicationDeterministic,
            ) 
        }
        
        if ( // host depletion is usually bulk of background
            params.pathogenDetection
                .qualityControl
                .hostDepletion 
        ) {
            reads = HostDepletion(
                reads, 
                databases.hostDepletion, 
                qualityControlParams.hostDepletionAligner
            ).out.reads
        }

        if ( // internal controls for coverage before further depletions
            params.pathogenDetection
                .qualityControl
                .internalControls 
        ) {
            reads = InternalControls(
                reads,
                databases.internalControls,
                qualityControlParams.internalControlsAligner
            ).out.reads
        }
        
        if ( // further background depletions of unwanted sequences
            params.pathogenDetection
                .qualityControl
                .backgroundDepletion
        ) {
            reads = BackgroundDepletion(
                reads, 
                databases.backgroundDepletion, 
                qualityControlParams.backgroundDepletionAligner
            ).out.reads
        }
        
        results = ReadScan.out.results.mix(
            SyntheticControls.out.results ?? Channel.empty(), 
            ReadQuality.out.results ?? Channel.empty(), 
            Deduplication.out.results  ?? Channel.empty(), 
            HostDepletion.out.results  ?? Channel.empty(), 
            InternalControls.out.results ?? Channel.empty(),
            BackgroundDepletion.out.results ?? Channel.empty()
        ) 
        | groupTuple

    emit:
        results = results

}