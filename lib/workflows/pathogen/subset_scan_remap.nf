include { ScrubbyReadsKrakenMinimapDepletionOnt as ScrubbyBackgroundDepletion } from '../../processes/scrubby' addParams(
    subdir: "viruses/depletion/background",
    result_file: "qc__scrubby__virus_background",
    kraken_taxa: params.virus_background_taxa,
    kraken_taxa_direct: params.virus_background_direct,
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
)
include { MinimapAlignPafOnt as MinimapAlign } from '../../processes/minimap2' addParams(
    subdir: "viruses/scan/alignments"
)
include { MinimapRealignBamOnt as MinimapRealign } from '../../processes/minimap2' addParams(
    subdir: "viruses/remap/alignments",
    covtobed_max_cov: params.taxa.alignment.viruses.remap.max_covtobed
)
include { VircovReferenceSelectionOnt as VircovReferenceSelection } from '../../processes/vircov' addParams(
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
include { VircovReferenceSelectionOnt as VircovReferenceSelectionBacteria } from '../../processes/vircov' addParams(
    subdir: "bacteria/scan",
    vircov_scan_min_mapq: params.taxa.alignment.bacteria.scan.min_mapq,
    vircov_scan_min_len: params.taxa.alignment.bacteria.scan.min_len,
    vircov_scan_min_cov: params.taxa.alignment.bacteria.scan.min_cov,
    vircov_scan_reads: params.taxa.alignment.bacteria.scan.min_reads,
    vircov_scan_coverage: params.taxa.alignment.bacteria.scan.min_coverage,
    vircov_scan_regions: params.taxa.alignment.bacteria.scan.min_regions,
    vircov_scan_regions_coverage: params.taxa.alignment.bacteria.scan.min_mapq,
    vircov_group_by: params.taxa.alignment.bacteria.scan.selection.group_by,
    vircov_group_sep: params.taxa.alignment.bacteria.scan.selection.group_sep,
    vircov_group_select_by: params.taxa.alignment.bacteria.scan.selection.select_by,
    vircov_group_select_segment_field: params.taxa.alignment.bacteria.scan.selection.segment_field,
    vircov_group_select_segment_field_nan: params.taxa.alignment.bacteria.scan.selection.segment_field_nan,
)
include { VircovReferenceSelectionOnt as VircovReferenceSelectionEukaryots } from '../../processes/vircov' addParams(
    subdir: "bacteria/scan",
    vircov_scan_min_mapq: params.taxa.alignment.eukaryots.scan.min_mapq,
    vircov_scan_min_len: params.taxa.alignment.eukaryots.scan.min_len,
    vircov_scan_min_cov: params.taxa.alignment.eukaryots.scan.min_cov,
    vircov_scan_reads: params.taxa.alignment.eukaryots.scan.min_reads,
    vircov_scan_coverage: params.taxa.alignment.eukaryots.scan.min_coverage,
    vircov_scan_regions: params.taxa.alignment.eukaryots.scan.min_regions,
    vircov_scan_regions_coverage: params.taxa.alignment.eukaryots.scan.min_mapq,
    vircov_group_by: params.taxa.alignment.eukaryots.scan.selection.group_by,
    vircov_group_sep: params.taxa.alignment.eukaryots.scan.selection.group_sep,
    vircov_group_select_by: params.taxa.alignment.eukaryots.scan.selection.select_by,
    vircov_group_select_segment_field: params.taxa.alignment.eukaryots.scan.selection.segment_field,
    vircov_group_select_segment_field_nan: params.taxa.alignment.eukaryots.scan.selection.segment_field_nan,
)
include { VircovRealignOnt as VircovRealign } from '../../processes/vircov' addParams(
    subdir: "viruses/remap",
    vircov_remap_min_mapq: params.taxa.alignment.viruses.remap.min_mapq,
    vircov_remap_min_len: params.taxa.alignment.viruses.remap.min_len,
    vircov_remap_min_cov: params.taxa.alignment.viruses.remap.min_cov,
    vircov_remap_reads: params.taxa.alignment.viruses.remap.min_reads,
    vircov_remap_coverage: params.taxa.alignment.viruses.remap.min_coverage,
    vircov_remap_regions: params.taxa.alignment.viruses.remap.min_regions,
    vircov_remap_regions_coverage: params.taxa.alignment.viruses.remap.min_mapq,
)
include { ConcatenateVircov } from '../../processes/vircov' addParams(
    subdir: "viruses/remap"
)
include { ScrubbyExtractVircovReadsOnt as ScrubbyExtractVircovReads } from '../../processes/scrubby' addParams(
    subdir: "viruses/remap/reads"
)
include { IvarConsensus } from '../../processes/ivar' addParams(
    subdir: "viruses/assembly/consensus"
)
include { ProcessVirusAlignment as ProcessAlignment } from '../../processes/cerebro' addParams(
    subdir: "viruses/"
)


/* 
==============================
DATABASE SUBSETS FOR ALIGNMENT
==============================
*/

include { MashScreenWinnerOnt as MashScreen } from '../../processes/mash' addParams(
    subdir: "alignment/mash_subset"
)
include { MashDatabaseSubset as MashDatabaseSubset } from '../../processes/cerebro' addParams(
    subdir: "alignment/mash_subset",
    min_shared_hashes: params.subset_min_shared_hashes
)
include { MinimapIndexSubsetOnt as MinimapIndexSubset } from '../../processes/minimap2' addParams(
    subdir: "alignment/mash_subset"
)

workflow background_depletion {
    take: 
        reads
        references                                   
        kraken_dbs                                                       
    main: 
        ScrubbyBackgroundDepletion(reads, kraken_dbs, references)       
    emit: 
        reads = ScrubbyBackgroundDepletion.out.reads
        results = ScrubbyBackgroundDepletion.out.results                 
}

workflow subset {
    take:
        reads                                                                                     
        mash_index                                                                    
        db_fasta
    main:
        MashScreen(reads, mash_index)                              
        MashDatabaseSubset(MashScreen.out, db_fasta)       
        MinimapIndexSubset(MashDatabaseSubset.out)          
    emit:
        reads = reads
        index = MinimapIndexSubset.out.index
        fasta = MinimapIndexSubset.out.fasta
}

workflow scan {
    take: 
        reads                                                 
        db_index
        db_fasta
        blacklist   
        domain                                 
    main: 
        // Align reads against viral database and select optimal references
        MinimapAlign(reads, db_index)   
        VircovReferenceSelection(MinimapAlign.out, db_fasta, blacklist)
        // For each sample, emit the selected references with their input reads individually for parallelisation
        // If there is just one output from the dynamic reference file outputs, it needs to be put into an
        // array - this does not automatically happen - see the additional mapping step
        aligned_references = VircovReferenceSelection.out.references | 
            map { data -> data[3] instanceof Collection ? data : tuple(data[0], data[1], data[2], [data[3]]) } |
            map { it[3].collect { el -> tuple(it[0], it[1], it[2], el) } } |
            flatMap
    emit:
        references = aligned_references // individual alignments and reference selections for paralellisation of realignments
        results = VircovReferenceSelection.out.results
}

workflow remap {
    take: 
        aligned_references      
        domain                                
    main: 
        // For each selected reference, realign the input reads and filter by alignment and coverage metrics
        MinimapRealign(aligned_references) 
        VircovRealign(MinimapRealign.out.reads, MinimapRealign.out.aligned)
        // Group outputs by sample and database name to concatenate results
        VircovRealign.out.aligned | groupTuple(by: [0, 1]) | ConcatenateVircov
        // ... and extract aligned reads for each reference 
        VircovRealign.out | ScrubbyExtractVircovReads
    emit:
        aligned = ScrubbyExtractVircovReads.out.aligned // id, db_name, ref_name, ref_fasta, bam, vircov, fwd_extracted, rev_extracted
        coverage = MinimapRealign.out.coverage          // id, db_name, covtxt
        results = ConcatenateVircov.out.results
}


workflow virus_consensus_assembly {
    take:
        aligned_data
    main:
        // Create a consensus reference sequence from aligned reads with iVar
        IvarConsensus(aligned_data)
    emit:
        IvarConsensus.out.consensus
}

workflow alignment {
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
        background_depletion(reads, references, kraken_dbs)

        if (database_subset) {
            subset(background_depletion.out.reads, mash_index, db_fasta)
            scan(subset.out.reads, subset.out.index, subset.out.fasta, blacklist, domain)
            remap(scan.out.references, domain)
        } else {
            scan(background_depletion.out.reads, db_index, db_fasta, blacklist, domain)
            remap(scan.out.references, domain)
        }
        
        // Optional consensus assembly
        if (consensus_assembly) { 
            consensus = virus_consensus_assembly(remap.out.aligned)
        } else {
            consensus = Channel.empty()
        }

        // Scan-remap-consensus result aggregation
        result_tables = scan.out.results.mix(remap.out.results) | groupTuple(by: [0, 1])

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
        ProcessAlignment(result_data)

    emit:
        results = ProcessAlignment.out.results
}   
