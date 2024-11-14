include { MinimapAlignPaf as MinimapAlignmentIllumina } from '../../processes/minimap2' addParams(
    subdir: "quality_control/phage/alignments"
)
include { ScrubbyAlignmentDepletion as AlignmentDepletionIllumina } from '../../processes/scrubby' addParams(
    subdir: "quality_control/phage/fastq",
    result_file: "qc__scrubby__phage",
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
)
include { VircovZero as AlignmentEvaluationIllumina } from '../../processes/vircov' addParams(
    subdir: "quality_control/phage",
    result_file: "qc__vircov__phage"
)

include { MinimapAlignPafOnt as MinimapAlignmentOnt } from '../../processes/minimap2' addParams(
    subdir: "quality_control/phage/alignments"
)
include { ScrubbyAlignmentDepletionOnt as AlignmentDepletionOnt } from '../../processes/scrubby' addParams(
    subdir: "quality_control/phage/fastq",
    result_file: "qc__scrubby__phage",
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
)
include { VircovZeroOnt as AlignmentEvaluationOnt } from '../../processes/vircov' addParams(
    subdir: "quality_control/phage",
    result_file: "qc__vircov__phage"
)

workflow phage_control {
    take: 
        reads                                         
        reference
        ont
    main: 
        if (ont) {
            MinimapAlignmentOnt(reads, reference)
            AlignmentDepletionOnt(MinimapAlignmentOnt.out)
            evaluation = AlignmentEvaluationOnt(MinimapAlignmentOnt.out, reference)
        } else {
            MinimapAlignmentIllumina(reads, reference)
            depletion = AlignmentDepletionIllumina(MinimapAlignmentIllumina.out)
            evaluation = AlignmentEvaluationIllumina(MinimapAlignmentIllumina.out, reference)
        }
    emit: 
        reads = depletion.reads
        results = depletion.results.mix(evaluation.results)                                     
}
