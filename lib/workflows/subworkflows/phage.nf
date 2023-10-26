include { MinimapAlignPAF as MinimapAlignment } from '../../processes/minimap2' addParams(
    subdir: "quality_control/phage/alignments"
)
include { ScrubbyAlignmentDepletion as AlignmentDepletion } from '../../processes/scrubby' addParams(
    subdir: "quality_control/phage/fastq",
    result_file: "qc__scrubby__phage",
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
)
include { VircovZero as AlignmentEvaluation } from '../../processes/vircov' addParams(
    subdir: "quality_control/phage",
    result_file: "qc__vircov__phage"
)

workflow phage_control {
    take: 
        reads                                         
        reference
    main: 
        MinimapAlignment(reads, reference)
        AlignmentDepletion(MinimapAlignment.out)
        AlignmentEvaluation(MinimapAlignment.out, reference)
    emit: 
        reads = AlignmentDepletion.out.reads
        results = AlignmentDepletion.out.results.mix(AlignmentEvaluation.out.results)                                     
}