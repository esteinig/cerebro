

include { from_reference_alignment_sample_sheet } from '../../utils';

include { MinimapAlignReferencePAF } from '../../processes/minimap2' addParams(
    subdir: "validation/reference_alignment"
)
include { VircovAlignReferenceZero } from '../../processes/vircov' addParams(
    subdir: "validation/reference_alignment"
)
include { FastpScan } from '../../processes/fastp' addParams(
    subdir: "validation/reference_alignment"
)


workflow reference_alignment_illumina {

    /* Workflow for simple alignment against sample specific reference genomes (Illumina) */
    
    take: 
        reads  // id, fwd, rev, fasta
    main: 
        reads | map { tuple(it[0], it[1], it[2]) } | FastpScan

        MinimapAlignReferencePAF(reads) | VircovAlignReferenceZero
    emit:
        VircovAlignReferenceZero.out.results

}

workflow reference_alignment {
    reads = from_reference_alignment_sample_sheet(params.validation.sample_sheet)
    reference_alignment_illumina(reads)
}