
include { quality_control_illumina }        from './subworkflows/quality_control';
include { quality_control_ont }             from './subworkflows/quality_control';
include { metagenome_assembly }             from './pathogen/metagenome_assembly';
include { kmer_profiling }                  from './pathogen/kmer_profiling';
include { alignment }                       from './pathogen/alignment';
include { alignment as panviral_alignment } from './pathogen/alignment';


workflow pathogen_detection {
    take:
        reads
        inputs
        ont
        pathogen
        panviral
    main:

        // ===========================
        // Quality control subworkflow
        // ===========================

        if (ont) {
            qc = quality_control_ont(
                reads, 
                inputs.ercc_fasta, 
                inputs.phage_fasta, 
                inputs.host_depletion_dbs, 
                inputs.host_depletion_references, 
                inputs.background_depletion_dbs, 
                inputs.background_depletion_references, 
                params.qc.host.depletion.enabled,
                params.qc.controls.phage.enabled,
                params.qc.background.depletion.enabled
            )
            qc_reads = qc[0]
        } else {
            qc = quality_control_illumina(
                reads, 
                inputs.adapter_fasta, 
                inputs.ercc_fasta, 
                inputs.phage_fasta, 
                inputs.host_depletion_dbs, 
                inputs.host_depletion_references, 
                inputs.background_depletion_dbs, 
                inputs.background_depletion_references, 
                params.qc.host.depletion.enabled,
                params.qc.controls.phage.enabled,
                params.qc.background.depletion.enabled,
                params.qc.deduplication.enabled,
                params.qc.deduplication.method
            )
            qc_reads = qc[0]
        }
       
        // ================================
        // Taxonomic classification modules
        // ================================

        if (pathogen) {

            log.info "Pathogen identification worklow enabled"

            if (params.taxa.alignment.enabled) {
                alignment_results = alignment(
                    qc_reads,
                    inputs.alignment_index,
                    inputs.alignment_fasta, 
                    null,
                    params.taxa.assembly.consensus.enabled,
                    false,
                    null,
                    ont
                )
            } else {
                alignment_results = Channel.empty()
            }

            if (params.taxa.kmer.enabeld) {
                kmer_profiling_results = kmer_profiling(
                    qc_reads,
                    inputs.kraken2_dbs,
                    inputs.metabuli_dbs,
                    ont
                )
            } else {
                kmer_profiling_results = Channel.empty()
            }
            
            qc_results = qc[1]

            result_files = qc_results.mix(
                alignment_results,
                kmer_profiling_results
            ) | groupTuple

        } else if (panviral) {

            log.info "Panviral enrichment subworkflow enabled"

            alignment_results = panviral_alignment(
                qc_reads,
                inputs.panviral_alignment_index,
                inputs.panviral_alignment_fasta, 
                inputs.panviral_blacklist,
                params.panviral.assembly.consensus.enabled,
                false,
                null,
                ont
            )
        }

    emit:
        results = result_files
}
