include { from_reference_alignment_sample_sheet } from '../../utils';

include { MinimapAlignReferencePAF } from '../../processes/minimap2' addParams(
    subdir: "validation/reference_alignment"
);
include { VircovAlignReferenceZero } from '../../processes/vircov' addParams(
    subdir: "validation/reference_alignment"
);
include { FastpScan } from '../../processes/fastp' addParams(
    subdir: "validation/reference_alignment"
);


workflow reference_alignment_illumina {

    /* Workflow for simple alignment against sample specific reference genomes */
    
    take: 
        reads   // id, fwd, rev, fasta
        inputs
    main: 
        if (params.validation.qc.enabled) {

            qc_reads = quality_control_illumina(
                reads | map { tuple(it[0], it[1], it[2]) }, 
                inputs.adapter_fasta, 
                inputs.ercc_fasta, 
                inputs.phage_fasta, 
                inputs.host_depletion_dbs, 
                inputs.host_depletion_references, 
                params.qc.deduplication.enabled,
                params.qc.deduplication.method,
                params.qc.host.depletion.enabled,
                params.qc.controls.phage.enabled
            ).out.reads

            align_reads = qc_reads.cross(reads) | map { tuple(it[0][0], it[0][1], it[0][2], it[1][3]) }
        
        } else {
            reads | map { tuple(it[0], it[1], it[2]) } | FastpScan
            align_reads = reads
        }

        MinimapAlignReferencePAF(align_reads) | VircovAlignReferenceZero
    emit:
        VircovAlignReferenceZero.out.results

}

workflow reference_alignment {
    
    take:
        inputs 
    main:
        reads = from_reference_alignment_sample_sheet(params.validation.sample_sheet)
        results = reference_alignment_illumina(reads, inputs)
    emit:
        results
}