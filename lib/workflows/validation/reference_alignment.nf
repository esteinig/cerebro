include { from_sample_sheet_reference_alignment_illumina; from_sample_sheet_reference_alignment_ont } from '../../utils';
include { quality_control_illumina; quality_control_ont } from '../subworkflows/quality_control';

include { MinimapAlignReferencePaf } from '../../processes/minimap2' addParams(
    subdir: "validation/reference_alignment"
);
include { MinimapAlignReferencePafOnt } from '../../processes/minimap2' addParams(
    subdir: "validation/reference_alignment"
);
include { VircovAlignReferenceZero } from '../../processes/vircov' addParams(
    subdir: "validation/reference_alignment"
);
include { VircovAlignReferenceZeroOnt } from '../../processes/vircov' addParams(
    subdir: "validation/reference_alignment"
);
include { FastpScan } from '../../processes/fastp' addParams(
    subdir: "validation/reference_alignment"
);
include { NanoqScan } from '../../processes/nanoq' addParams(
    subdir: "validation/reference_alignment",
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
            ).reads

            align_reads = qc_reads.cross(reads) | map { tuple(it[0][0], it[0][1], it[0][2], it[1][3]) }
        
        } else {
            reads | map { tuple(it[0], it[1], it[2]) } | FastpScan
            align_reads = reads
        }

        MinimapAlignReferencePaf(align_reads) | VircovAlignReferenceZero
    emit:
        VircovAlignReferenceZero.out.results

}


workflow reference_alignment_ont {

    /* Workflow for simple alignment against sample specific reference genomes */
    
    take: 
        reads   // id, fastq, fasta
        inputs
    main: 
        if (params.validation.qc.enabled) {

            qc_reads = quality_control_ont(
                reads | map { tuple(it[0], it[1]) }, 
                inputs.adapter_fasta, 
                inputs.ercc_fasta, 
                inputs.phage_fasta, 
                inputs.host_depletion_dbs, 
                inputs.host_depletion_references, 
                params.qc.deduplication.enabled,
                params.qc.deduplication.method,
                params.qc.host.depletion.enabled,
                params.qc.controls.phage.enabled
            ).reads

            align_reads = qc_reads.cross(reads) | map { tuple(it[0][0], it[0][1], it[1][3]) }
        
        } else {
            reads | map { tuple(it[0], it[1]) } | NanoqScan
            align_reads = reads
        }

        MinimapAlignReferencePaf(align_reads) | VircovAlignReferenceZero
    emit:
        VircovAlignReferenceZero.out.results

}

workflow reference_alignment {
    
    take:
        inputs 
        ont
    main:
        if (ont) {
            reads = from_sample_sheet_reference_alignment_illumina(params.validation.sample_sheet)
            results = reference_alignment_illumina(reads, inputs)
        } else {
            reads = from_sample_sheet_reference_alignment_ont(params.validation.sample_sheet)
            results = reference_alignment_ont(reads, inputs)
        }
    emit:
        results
}