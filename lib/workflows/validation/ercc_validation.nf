/*
    ================================================================================================
    ERCC/EDCC de-duplication and biomass calculation experiment for assay validation [Validation #1]
    ================================================================================================


    Validate the use of ERCC and our EDCC protocols for accurate biomass calculations of host 
    and total biomass inputs, particularly in later experiments on the limit of detection for 
    various organisms.

    =================
    ** Description **
    =================

    We assume that the input is a series of known human cell biomass libraries (100 pg, 75 pg, 50 pg ...)
    with different ERCC/EDCC input mass (25 pg, 15 pg, 10 pg, 5 pg ...) using the NEB Dual Index UMI
    protocol as described in the assay protocol, ideally with technical replicates.

    Since the quantification of biomass requires deduplication of ERCC and other reads, we 
    test various approaches for the most useful implementation that allows us to accurately
    back-calculate the human cell biomass of the input libraries. 

    Deduplication with UMIs is necessary because any naive deduplication of ERCCs
    would deplete the initially duplicated input sequences themselves, leading to
    overestimates of biomass.

    ==============
    ** Workflow **
    ==============

    In general the quality control pipeline runs the following steps:


        ==> UMI Deduplication -> ERCC/EDCC alignment and depletion -> Read QC -> Human alignment and depletion


    implemented in:

        - `deduplication_naive_umi`
            
          UMI deduplication using naive clustering of forward sequences with UMI prepended by exact sequence identity

        - `deduplication_calib`

          Reference free deduplication using Calib - slightly more involved computations as outlined in their paper

        - `deduplication_umi_tools`
        
          Reference alignment and position based deduplication with UMI-tools (host, ercc) - standard for transcript 
          quanitification and single cell, but may not be suitable because of metagenomic databases and runtimes.

        - `deduplication_umi_tools_naive`
            
          Reference alignment and position based deduplication with umi-tools followed by naive deduplication of remaining reads.
          Combines the alignment based deduplication of single references () with naive reference free approach for remaining
          sequence content (i.e. metagenomic content, microbial taxa)


    Variations of this are included as controls to assess expected over- or underestimates of host biomass:


        ==> ERCC/EDCC alignment and depletion -> Read QC -> Human alignment and depletion


        - `deduplication_none`
        
          Does not run deduplication and therefore retains effect of PCR amplification bias between ERCC and human input.
        

        ==> Deduplication -> ERCC/EDCC alignment and depletion -> Read QC -> Human alignment and depletion


        - `deduplication_naive_before`

        
          Runs naive deduplication before ercc removal: we expect overstimate of host biomass because original input ercc sequences are deduplicated
        

        ==> ERCC/EDCC alignment and depletion -> Deduplication -> Read QC -> Human alignment and depletion


        - `deduplication_naive_after`

          Runs naive deduplication after ercc removal: we expect underestimate of host biomass because ercc is not deduplicated but host sequences are 


    Note that reference based deduplication with UMI-tools takes a logn time to compute due to the alignment of
    host sequences! It may not be suitable for production implementation.
    

    ============
    ** Tools **
    ============
    
    Tools used in this workflow are:
        
        - cerebro dedup-naive 
        - calib
        - umi_tools
        - minimap2
        - fastp
        - scrubby
          - kraken2
          - strobealign
          - bowtie2
          - bwa

    Input for the UMI-tools step must be an index of the host reference sequences that includes the ERCC (`params.host_ercc_index`)

*/


include { quality_control as dedup_none  }            from "../quality_control";
include { quality_control as dedup_fastp }            from "../quality_control";
include { quality_control as dedup_naive }            from "../quality_control";
include { quality_control as dedup_umi_naive }        from "../quality_control";
include { quality_control as dedup_umi_calib }        from "../quality_control";
include { quality_control as dedup_umi_tools }        from "../quality_control";
include { quality_control as dedup_umi_tools_naive }  from "../quality_control";


workflow ercc_validation {

    take: 
        reads
        adapter_fasta
        ercc_fasta
        kraken_dbs
        host_references
        host_ercc_index
    main: 

        // Note that this may generate a large amount of intermediate file space
        // we are waiting for the cleanup implementation in the next major 
        // Nextflow release to be able to run these experiments in a single
        
        if (params.validation.ercc.deduplication.none) {
          dedup_none(reads, adapter_fasta, ercc_fasta, kraken_dbs, host_references, host_ercc_index, "none")
        }
        if (params.validation.ercc.deduplication.fastp) {
          dedup_none(reads, adapter_fasta, ercc_fasta, kraken_dbs, host_references, host_ercc_index, "fastp")
        }
        if (params.validation.ercc.deduplication.naive) {
          dedup_naive_before(reads, adapter_fasta, ercc_fasta, kraken_dbs, host_references, host_ercc_index, "naive")
        }
        if (params.validation.ercc.deduplication.umi_naive) {
          dedup_umi_naive(reads, adapter_fasta, ercc_fasta, kraken_dbs, host_references, host_ercc_index, "umi-naive")
        }
        if (params.validation.ercc.deduplication.umi_calib) {
          dedup_umi_calib(reads, adapter_fasta, ercc_fasta, kraken_dbs, host_references, host_ercc_index, "umi-calib")
        }
        if (params.validation.ercc.deduplication.umi_tools) {
          dedup_umi_tools(reads, adapter_fasta, ercc_fasta, kraken_dbs, host_references, host_ercc_index, "umi-tools")
        }
        if (params.validation.ercc.deduplication.umi_tools_naive) {
          dedup_umi_tools_naive(reads, adapter_fasta, ercc_fasta, kraken_dbs, host_references, host_ercc_index, "umi-tools-naive")
        }

}