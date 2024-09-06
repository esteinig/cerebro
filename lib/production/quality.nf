/* ==============================
 * Quality control for production
 * ==============================
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
} from "../processes/quality"

workflow QualityControl {
    take:
        reads
        databases
    main:

        InputScan(reads)

        // Quality control background only or start with synthetic controls
        if (params.cerebroConfig.qualityControlBackgroundOnly) {
            reads = BackgroundDepletion(
                reads, 
                databases.backgroundDepletion, 
                params.qualityControl.backgroundDepletionAligner
            ).reads

        } else {
            if (params.qualityControl.syntheticControls) {
                reads = SyntheticControls(
                    reads,
                    databases.syntheticControls,
                    params.qualityControl.syntheticControlsAligner
                ).reads
            }
        }

        // Deduplicate before read quality control or not
        if (params.cerebroConfig.qualityControlDeduplicateFirst) {

            if (params.qualityControl.readDeduplication) {
                reads = Deduplication(
                    reads,
                    params.qualityControl.readDeduplicationHead,
                    params.qualityControl.readDeduplicationDeterministic
                ).reads
            }
            if (params.qualityControl.readQuality) {
                reads = ReadQuality(reads).reads
            }
        } else {
            if (params.qualityControl.readQuality) {
                reads = ReadQuality(reads).reads
            }
            if (params.qualityControl.readDeduplication) {
                reads = Deduplication(
                    reads,
                    params.qualityControl.readDeduplicationHead,
                    params.qualityControl.readDeduplicationDeterministic
                ).reads
            }
        }


        // Complete optional alignment depletions
        if (!params.cerebroConfig.qualityControlBackgroundOnly) {
            if (params.qualityControl.hostDepletion) {
                reads = HostDepletion(
                    reads, 
                    databases.hostDepletion, 
                    params.qualityControl.hostDepletionAligner
                ).reads
            }

            if (params.qualityControl.internalControls) {
                reads = InternalControls(
                    reads,
                    databases.internalControls,
                    params.qualityControl.internalControlsAligner
                ).reads
            }

            if (params.qualityControl.backgroundDepletion) {
                reads = BackgroundDepletion(
                    reads, 
                    databases.backgroundDepletion, 
                    params.qualityControl.backgroundDepletionAligner
                ).reads
            }
        }

        OutputScan(reads)
        
        // collect results by sample id
        results = InputScan.out.results.mix(
            params.qualityControl.syntheticControls   ? SyntheticControls.out.results   : Channel.empty(), 
            params.qualityControl.readQuality         ? ReadQuality.out.results         : Channel.empty(), 
            params.qualityControl.readDeduplication   ? Deduplication.out.results       : Channel.empty(), 
            params.qualityControl.hostDepletion       ? HostDepletion.out.results       : Channel.empty(), 
            params.qualityControl.internalControl     ? InternalControls.out.results    : Channel.empty(),
            params.qualityControl.backgroundDepletion ? BackgroundDepletion.out.results : Channel.empty(),
            OutputScan.out.results
        )

        // prcoess results to json and get qc tables
        json    = results | groupTuple | ProcessOutput 
        tables  = json    | collect    | QualityControlTables

    emit:
        tables  = tables
        json    = json
        results = results
        reads   = reads
}
