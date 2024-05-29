

include { MinimapAlignReferencePAF } from '../processes/minimap2';
include { VircovAlignReferenceZero } from '../processes/vircov';
include { from_reference_alignment_sample_sheet } from '../../utils';



workflow reference_alignment_illumina {

    /* Workflow for simple alignment against sample specific reference genomes (Illumina) */
    
    take: 
        reads  // id, fwd, rev, fasta
    main: 
        MinimapAlignReferencePAF(reads) | VircovAlignReferenceZero
    emit:
        VircovAlignReferenceZero.out.results

}

workflow {
    reads = from_reference_alignment_sample_sheet(params.sample_sheet)
    reference_alignment_illumina(reads)
}