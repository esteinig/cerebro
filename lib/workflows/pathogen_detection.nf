
include { eukaryots_subset_alignment } from './pathogen/subset_alignment'
include { bacteria_subset_alignment }  from './pathogen/subset_alignment'
include { virus_detection }            from './pathogen/virus_detection'
include { kmer_profiling }             from './pathogen/kmer_profiling'
include { meta_assembly }              from './pathogen/meta_assembly'
include { quality_control }            from './subworkflows/quality_control'

workflow pathogen_detection {
    take:
        reads
        inputs
    main:
        // ===========================
        // Quality control subworkflow
        // ===========================

        quality_control(
            reads, 
            inputs.adapter_fasta, 
            inputs.ercc_fasta, 
            inputs.phage_fasta, 
            inputs.host_depletion_dbs, 
            inputs.host_depletion_references, 
            inputs.host_ercc_index, 
            params.qc.deduplication.enabled,
            params.qc.deduplication.method,
            params.qc.host.depletion.enabled,
            params.qc.controls.phage.enabled
        )
        
        // ================================
        // Taxonomic classification modules
        // ================================

        // =================
        // Alignment modules
        // =================

        // Viral profiling, coverage assessment, consensus assembly after background depletion
        if (params.taxa.alignment.enabled && params.taxa.alignment.viruses.enabled) {
            virus_detection(
                quality_control.out.reads, 
                inputs.virus_background_references, 
                inputs.virus_background_dbs, 
                inputs.virus_db_index, 
                inputs.virus_db_fasta, 
                inputs.virus_blacklist,
                params.taxa.assembly.viruses.enabled
            )
            align_virus_results = virus_detection.out.results
        } else {
            align_virus_results = Channel.empty()
        }

        // Bacterial subset alignments
        if (params.taxa.alignment.enabled && params.taxa.alignment.bacteria.enabled) {
            bacteria_subset_alignment(
                quality_control.out.reads, 
                inputs.bacteria_mash_index, 
                inputs.bacteria_fasta
            )
            align_bacteria_results = bacteria_subset_alignment.out.results
        } else {
            align_bacteria_results = Channel.empty()
        }

        // Eukaryotic subset alignments
        if (params.taxa.alignment.enabled && params.taxa.alignment.eukaryots.enabled){
            eukaryots_subset_alignment(
                quality_control.out.reads, 
                inputs.eukaryots_mash_index, 
                inputs.eukaryots_fasta
            )
            align_eukaryots_results = eukaryots_subset_alignment.out.results
        } else {
            align_eukaryots_results = Channel.empty()
        }
        
        // =============
        // K-mer modules
        // =============

        if (params.taxa.kmer.enabled && params.taxa.kmer.kraken2.enabled){
            kmer_profiling(
                quality_control.out.reads, 
                inputs.kraken2_dbs
            )
            kmer_results = kmer_profiling.out.results
        } else {
            kmer_results = Channel.empty()
        }
        
        // ================
        // Assembly modules
        // ================

        if (params.taxa.assembly.enabled && params.taxa.assembly.meta.enabled) {
            meta_assembly(
                quality_control.out.reads, 
                inputs.meta_blast_nt, 
                inputs.meta_diamond_nr
            )
            meta_assembly_results = meta_assembly.out.results
        } else {
            meta_assembly_results = Channel.empty()
        }

        // We wait until all processes have finished by collecting all their results file outputs
        result_files = quality_control.out.results.mix(
            align_virus_results, 
            align_bacteria_results, 
            align_eukaryots_results, 
            kmer_results, 
            meta_assembly_results
        ) | groupTuple

    emit:
        results = result_files
}
           