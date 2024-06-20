
include { quality_control_illumina }      from './subworkflows/quality_control';
include { quality_control_ont }           from './subworkflows/quality_control';
include { metagenome_assembly }           from './pathogen/metagenome_assembly';
include { kmer_profiling }                from './pathogen/kmer_profiling';
include { alignment }                     from './pathogen/alignment';


workflow pathogen_detection {
    take:
        reads
        inputs
        ont
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


        // +++++++++++++++++++++
        // Alignment subworkflow
        // +++++++++++++++++++++

        if (params.taxa.alignment.enabled) {
            alignment_results = alignment(
                qc_reads,
                inputs.virus_background_references,
                inputs.virus_background_dbs,
                inputs.virus_db_index,
                inputs.virus_db_fasta, 
                inputs.virus_blacklist,
                params.taxa.assembly.consensus.enabled,
                false,
                null,
                "viruses",
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

    emit:
        results = result_files
}



        // // Bacterial subset alignments
        // if (params.taxa.alignment.enabled && params.taxa.alignment.domains.bacteria.enabled) {
        //     bacteria_subset_alignment(
        //         qc_reads, 
        //         inputs.bacteria_mash_index, 
        //         inputs.bacteria_fasta
        //     )
        //     align_bacteria_results = bacteria_subset_alignment.out.results
        // } else {
        //     align_bacteria_results = Channel.empty()
        // }

        // // Eukaryotic subset alignments
        // if (params.taxa.alignment.enabled && params.taxa.alignment.domains.eukaryots.enabled){
        //     eukaryots_subset_alignment(
        //         qc_reads, 
        //         inputs.eukaryots_mash_index, 
        //         inputs.eukaryots_fasta
        //     )
        //     align_eukaryots_results = eukaryots_subset_alignment.out.results
        // } else {
        //     align_eukaryots_results = Channel.empty()
        // }
        
        // =============
        // K-mer modules
        // =============

        // if (params.taxa.kmer.enabled && params.taxa.kmer.kraken2.enabled && !params.taxa.kmer.bracken.enabled){
        //     // Pathogen detection as per Nature Protocols (Kraken2Uniq)
        //     kmer_pathogen_detection(
        //         qc_reads, 
        //         inputs.kraken2_dbs,
        //         ont
        //     )
        //     kmer_results = kmer_pathogen_detection.out.results
        // } else if (params.taxa.kmer.enabled && params.taxa.kmer.kraken2.enabled && params.taxa.kmer.bracken.enabled) {
        //     // Pathogen profiling with Kraken2 + Bracken
        //     kmer_pathogen_profiling(
        //         qc_reads, 
        //         inputs.kraken2_dbs,
        //         ont
        //     )
        //     kmer_results = kmer_pathogen_profiling.out.results
        // } else {
        //     kmer_results = Channel.empty()
        // }
        
        // // ================
        // // Assembly modules
        // // ================

        // if (params.taxa.assembly.enabled && params.taxa.assembly.meta.enabled) {
        //     metagenome_assembly(
        //         qc_reads, 
        //         inputs.meta_blast_nt, 
        //         inputs.meta_diamond_nr
        //     )
        //     meta_assembly_results = metagenome_assembly.out.results
        // } else {
        //     meta_assembly_results = Channel.empty()
        // }


        // // We wait until all processes have finished by collecting all their results file outputs
        // result_files = qc_results.mix(
        //     align_virus_results, 
        //     align_bacteria_results, 
        //     align_eukaryots_results, 
        //     kmer_results, 
        //     meta_assembly_results
        // ) | groupTuple