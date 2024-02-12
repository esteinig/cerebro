
include { kmer_pathogen_profiling as kmer_pathogen_profiling_ont }    from './pathogen/kmer_profiling'
include { metagenome_assembly }        from './pathogen/metagenome_assembly'
include { quality_control_ont }        from './subworkflows/quality_control'
include { quality_control_illumina }   from './subworkflows/quality_control'

workflow culture_identification {
    take:
        ont_reads
        illumina_reads
        inputs
    main:

        // ===========================
        // Quality control subworkflow
        // ===========================

        qc_ont = quality_control_ont(
            ont_reads, 
            inputs.ercc_fasta, 
            inputs.phage_fasta, 
            inputs.host_depletion_dbs, 
            inputs.host_depletion_references, 
            params.qc.host.depletion.enabled,
            params.qc.controls.phage.enabled
        )

        qc_illumina = quality_control_illumina(
            illumina_reads, 
            inputs.adapter_fasta, 
            inputs.ercc_fasta, 
            inputs.phage_fasta, 
            inputs.host_depletion_dbs, 
            inputs.host_depletion_references, 
            params.qc.deduplication.enabled,
            params.qc.deduplication.method,
            params.qc.host.depletion.enabled,
            params.qc.controls.phage.enabled
        )

        // ====================
        // K-mer read profiling
        // ====================
        
        // Separately for Illumina and ONT
        if (params.taxa.kmer.enabled && params.taxa.kmer.kraken2.enabled && params.taxa.kmer.bracken.enabled) {
            taxa_ont = kmer_pathogen_profiling_ont(qc_ont.reads, inputs.kraken2_dbs, true)
            taxa_illumina = kmer_pathogen_profiling_illumina(qc_illumina.reads, inputs.kraken2_dbs, false)

            kmer_results = taxa_ont.results.mix(taxa_illumina.results)
        } else {
            kmer_results = Channel.empty()
        }

        // ============================
        // Hybrid assembly - Dragonflye
        // ============================

        qc_ont.reads.mix(qc_illumina.reads) | groupTuple | view
        
}