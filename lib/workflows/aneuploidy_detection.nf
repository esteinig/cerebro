include { MinimapAneuploidy } from '../processes/minimap2' addParams(
    subdir: 'host/aneuploidy/alignments'
)
include { CnvKitAneuploidy } from '../processes/cnvkit' addParams(
    subdir: 'host/aneuploidy/copy_number_variation'
)

include { quality_control } from './subworkflows/quality_control';

// Aneuploidy detection through whole genome copy number variation
workflow aneuploidy_detection {
    take:
        reads
        inputs
    main:
        // Quality control of libraries without host removal
        quality_control(
            reads, 
            inputs.adapter_fasta, 
            inputs.ercc_fasta, 
            inputs.phage_fasta, 
            inputs.host_depletion_dbs, 
            inputs.host_depletion_references, 
            inputs.host_ercc_index, 
            false,                                      // primary deduplication disabled, uses host alignment deduplication instead
            null,
            false,                                      // host removal disabled 
            params.qc.controls.phage.enabled
        )

        aneuploidy_cnv(
            quality_control.out.reads, 
            inputs
        )

    emit:
        results = aneuploidy_cnv.out.results
}

workflow aneuploidy_cnv {
    take:
        reads
        inputs
    main:
        // Align reads against human reference [sorted BAM]
        MinimapAneuploidy(reads, inputs.aneuploidy_reference_index)

        // Use precomputed BAM control alignment against same reference
        // to run copy number variation detection from whole genome
        if (params.host.aneuploidy.cnvkit.enabled) {
            CnvKitAneuploidy(MinimapAneuploidy.out.alignment, inputs.aneuploidy_controls, inputs.aneuploidy_reference_fasta)
            results = CnvKitAneuploidy.out.results | groupTuple
        } else {
            results = Channel.empty()
        }
    emit: 
        results = results
}
