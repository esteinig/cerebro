include { ScrubbyReadsKrakenMinimapDepletionOnt as ScrubbyBackgroundDepletionONT } from '../../processes/scrubby' addParams(
    result_file: "qc__scrubby__virus_background",
    kraken_taxa: params.virus_background_taxa,
    kraken_taxa_direct: params.virus_background_direct,
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
)
include { MinimapAlignPafOnt as MinimapAlignONT } from '../../processes/minimap2' addParams(
    subdir: "viruses/scan/alignments"
)
include { MinimapRealignBamOnt as MinimapRealignONT } from '../../processes/minimap2' addParams(
    subdir: "viruses/remap/alignments",
    covtobed_max_cov: params.taxa.alignment.viruses.remap.max_covtobed
)
include { VircovReferenceSelectionOnt as VircovReferenceSelectionONT } from '../../processes/vircov' addParams(
    subdir: "viruses/scan",
    vircov_scan_min_mapq: params.taxa.alignment.viruses.scan.min_mapq,
    vircov_scan_min_len: params.taxa.alignment.viruses.scan.min_len,
    vircov_scan_min_cov: params.taxa.alignment.viruses.scan.min_cov,
    vircov_scan_reads: params.taxa.alignment.viruses.scan.min_reads,
    vircov_scan_coverage: params.taxa.alignment.viruses.scan.min_coverage,
    vircov_scan_regions: params.taxa.alignment.viruses.scan.min_regions,
    vircov_scan_regions_coverage: params.taxa.alignment.viruses.scan.min_mapq,
    vircov_group_by: params.taxa.alignment.viruses.scan.selection.group_by,
    vircov_group_sep: params.taxa.alignment.viruses.scan.selection.group_sep,
    vircov_group_select_by: params.taxa.alignment.viruses.scan.selection.select_by,
    vircov_group_select_segment_field: params.taxa.alignment.viruses.scan.selection.segment_field,
    vircov_group_select_segment_field_nan: params.taxa.alignment.viruses.scan.selection.segment_field_nan,
)
include { VircovRealignOnt as VircovRealignONT } from '../../processes/vircov' addParams(
    subdir: "viruses/remap",
    vircov_remap_min_mapq: params.taxa.alignment.viruses.remap.min_mapq,
    vircov_remap_min_len: params.taxa.alignment.viruses.remap.min_len,
    vircov_remap_min_cov: params.taxa.alignment.viruses.remap.min_cov,
    vircov_remap_reads: params.taxa.alignment.viruses.remap.min_reads,
    vircov_remap_coverage: params.taxa.alignment.viruses.remap.min_coverage,
    vircov_remap_regions: params.taxa.alignment.viruses.remap.min_regions,
    vircov_remap_regions_coverage: params.taxa.alignment.viruses.remap.min_mapq,
)
include { ConcatenateVircov as ConcatenateVircovONT } from '../../processes/vircov' addParams(
    subdir: "viruses/remap"
)
include { ScrubbyExtractVircovReadsOnt as ScrubbyExtractVircovReadsONT } from '../../processes/scrubby' addParams(
    subdir: "viruses/remap/reads"
)
include { IvarConsensus as IvarConsensusONT } from '../../processes/ivar' addParams(
    subdir: "viruses/assembly/consensus"
)
include { ProcessVirusAlignment as ProcessAlignmentONT } from '../../processes/cerebro' addParams(
    subdir: "viruses/"
)


/* 
==============================
DATABASE SUBSETS FOR ALIGNMENT
==============================
*/

include { MashScreenWinnerOnt as MashScreenONT } from '../../processes/mash' addParams(
    subdir: "alignment/mash_subset"
)
include { MashDatabaseSubset as MashDatabaseSubsetONT } from '../../processes/cerebro' addParams(
    subdir: "alignment/mash_subset",
    min_shared_hashes: params.subset_min_shared_hashes
)
include { MinimapIndexSubsetOnt as MinimapIndexSubsetONT } from '../../processes/minimap2' addParams(
    subdir: "alignment/mash_subset"
)

workflow background_depletion_ont {
    take: 
        reads
        references                                   
        kraken_dbs    
        domain                                                   
    main: 
        ScrubbyBackgroundDepletionONT(reads, kraken_dbs, references, domain)       
    emit: 
        reads = ScrubbyBackgroundDepletionONT.out.reads
        results = ScrubbyBackgroundDepletionONT.out.results                 
}

workflow subset_ont {
    take:
        reads                                                                                     
        mash_index                                                                    
        db_fasta
    main:
        MashScreenONT(reads, mash_index)                              
        MashDatabaseSubsetONT(MashScreenONT.out, db_fasta)       
        MinimapIndexSubsetONT(MashDatabaseSubsetONT.out)          
    emit:
        reads = reads
        index = MinimapIndexSubsetONT.out.index
        fasta = MinimapIndexSubsetONT.out.fasta
}

workflow scan_ont {
    take: 
        reads                                                 
        db_index
        db_fasta
        blacklist   
        domain                                 
    main: 
        // Align reads against viral database and select optimal references
        MinimapAlignONT(reads, db_index)   
        VircovReferenceSelectionONT(MinimapAlignONT.out, db_fasta, blacklist)
        // For each sample, emit the selected references with their input reads individually for parallelisation
        // If there is just one output from the dynamic reference file outputs, it needs to be put into an
        // array - this does not automatically happen - see the additional mapping step
        aligned_references = VircovReferenceSelectionONT.out.references | 
            map { data -> data[3] instanceof Collection ? data : tuple(data[0], data[1], data[2], [data[3]]) } |
            map { it[3].collect { el -> tuple(it[0], it[1], it[2], el) } } |
            flatMap
    emit:
        references = aligned_references // individual alignments and reference selections for paralellisation of realignments
        results = VircovReferenceSelectionONT.out.results
}

workflow remap_ont {
    take: 
        aligned_references      
        domain                                
    main: 
        // For each selected reference, realign the input reads and filter by alignment and coverage metrics
        MinimapRealignONT(aligned_references) 
        VircovRealignONT(MinimapRealignONT.out.reads, MinimapRealignONT.out.aligned)
        // Group outputs by sample and database name to concatenate results
        VircovRealignONT.out.aligned | groupTuple(by: [0, 1]) | ConcatenateVircovONT
        // ... and extract aligned reads for each reference 
        VircovRealignONT.out | ScrubbyExtractVircovReadsONT
    emit:
        aligned = ScrubbyExtractVircovReadsONT.out.aligned // id, db_name, ref_name, ref_fasta, bam, vircov, fwd_extracted, rev_extracted
        coverage = MinimapRealignONT.out.coverage          // id, db_name, covtxt
        results = ConcatenateVircovONT.out.results
}


workflow virus_consensus_assembly_ont {
    take:
        aligned_data
    main:
        // Create a consensus reference sequence from aligned reads with iVar
        IvarConsensusONT(aligned_data)
    emit:
        IvarConsensusONT.out.consensus
}

workflow alignment_ont {
    take:
        reads
        references
        kraken_dbs
        db_index
        db_fasta
        blacklist
        consensus_assembly
        database_subset
        mash_index
        domain
    main:

        // Background depletion
        background_depletion_ont(reads, references, kraken_dbs)

        if (database_subset) {
            subset_ont(background_depletion_ont.out.reads, mash_index, db_fasta)
            scan_ont(subset_ont.out.reads, subset.out.index, subset.out.fasta, blacklist, domain)
            remap_ont(scan_ont.out.references, domain)
        } else {
            scan_ont(background_depletion_ont.out.reads, db_index, db_fasta, blacklist, domain)
            remap_ont(scan_ont.out.references, domain)
        }
        
        // Optional consensus assembly
        if (consensus_assembly) { 
            consensus = virus_consensus_assembly_ont(remap.out.aligned)
        } else {
            consensus = Channel.empty()
        }

        // Scan-remap-consensus result aggregation
        result_tables = scan_ont.out.results.mix(remap.out.results) | groupTuple(by: [0, 1])

        // Mixing outputs does not guarantee order, remap results are not guaranteed, 
        // we will account for this inside the process
        result_data = result_tables.mix(
            consensus | groupTuple(by: [0, 1]),
            remap.out.coverage | groupTuple(by: [0, 1])
        ) | groupTuple(by: [0, 1]) | map { 
            data -> tuple(
                data[0], 
                data[1], 
                data[2][0] instanceof Collection ? data[2][0] : [], 
                data[2][1] instanceof Collection ? data[2][1] : [],
                data[2][2] instanceof Collection ? data[2][2] : []
            )
        } // check for optional outputs and account for missing data (e.g. as in negative controls)

        // Aggregation of scan-remap-consensus results
        ProcessAlignmentONT(result_data)

    emit:
        results = ProcessAlignmentONT.out.results
}   
