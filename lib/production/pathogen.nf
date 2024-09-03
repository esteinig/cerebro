/* ==================================
 * PATHOGEN DETECTION FOR PRODUCTION
 * ==================================
*/

include { 
    InputScan
    SyntheticControls
    ReadQuality
    Deduplication
    HostDepletion
    InternalControls
    BackgroundDepletion
    OutputScan
    ProcessOutput
    QualityControlTables
} from "../processes/pathogen"


workflow PathogenDetection {

    take:
        reads
        qualityControlDatabases
    main:
        qualityControl = QualityControl(reads, qualityControlDatabases)

        qualityControlResults = qualityControl.results | groupTuple

        qualityControlJson = qualityControlResults | ProcessOutput 
        qualityControlJson | collect | QualityControlTables



    emit:
        qc = qualityControlResults

}

workflow QualityControl {
    take:
        reads
        databases
    main:

        qualityControlParams = params.pathogenDetection.qualityControl

        // total read counts and table before depletions
        InputScan(reads)

        // first for synthetic controls for accurate biomass estimates
        if (qualityControlParams.syntheticControls) {
            reads = SyntheticControls(
                reads,
                databases.syntheticControls,
                qualityControlParams.syntheticControlsAligner
            ).reads
        }

        // second for read quality control and trimming
        if (qualityControlParams.readQuality) {
            reads = ReadQuality(reads).reads
        }

        // third for read head + umi deduplication
        if (qualityControlParams.readDeduplication) {
            reads = Deduplication(
                reads,
                qualityControlParams.readDeduplicationHead,
                qualityControlParams.readDeduplicationDeterministic
            ).reads
        }
        
        // host depletion is usually bulk of background
        if (qualityControlParams.hostDepletion) {
            reads = HostDepletion(
                reads, 
                databases.hostDepletion, 
                qualityControlParams.hostDepletionAligner
            ).reads
        }

        // internal controls for coverage before further depletions
        if (qualityControlParams.internalControls) {
            reads = InternalControls(
                reads,
                databases.internalControls,
                qualityControlParams.internalControlsAligner
            ).reads
        }

        // further background depletions of unwanted sequences
        if (qualityControlParams.backgroundDepletion) {
            reads = BackgroundDepletion(
                reads, 
                databases.backgroundDepletion, 
                qualityControlParams.backgroundDepletionAligner
            ).reads
        }

        // total read counts and table after depletions
        OutputScan(reads)
        
        // collect results by sample id
        results = InputScan.out.results.mix(
            qualityControlParams.syntheticControls   ? SyntheticControls.out.results   : Channel.empty(), 
            qualityControlParams.readQuality         ? ReadQuality.out.results         : Channel.empty(), 
            qualityControlParams.readDeduplication   ? Deduplication.out.results       : Channel.empty(), 
            qualityControlParams.hostDepletion       ? HostDepletion.out.results       : Channel.empty(), 
            qualityControlParams.internalControl     ? InternalControls.out.results    : Channel.empty(),
            qualityControlParams.backgroundDepletion ? BackgroundDepletion.out.results : Channel.empty(),
            OutputScan.out.results
        )

    emit:
        results = results

}


// workflow TaxonomicProfile {
//     take:
//         reads
//         databases
//     main:
        

//     emit:

// }


// workflow MetagenomeAssembly {
//     take:
//         reads
//         databases
//     main:


//     emit:
// }