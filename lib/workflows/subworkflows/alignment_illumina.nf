// ILLUMINA

include { ScrubbyReadsKrakenMinimapDepletion as ScrubbyBackgroundDepletion } from '../../processes/scrubby' addParams(
    result_file: "qc__scrubby__virus_background",
    subdir: "pathogens/alignment/background_depletion",
    kraken_taxa: params.qc.background.depletion.taxa,
    kraken_taxa_direct: params.qc.background.depletion.direct,
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
);
include { MinimapAlignPaf as MinimapAlign } from '../../processes/minimap2' addParams(
    subdir: "pathogens/alignment/scan"
);
include { MinimapRealignBam as MinimapRealign } from '../../processes/minimap2' addParams(
    subdir: "pathogens/alignment/remap"
);
include { VircovReferenceSelection as VircovReferenceSelection } from '../../processes/vircov' addParams(
    subdir: "pathogens/alignment/scan",
    vircov_scan_min_mapq: params.taxa.alignment.scan.min_mapq,
    vircov_scan_min_len: params.taxa.alignment.scan.min_len,
    vircov_scan_min_cov: params.taxa.alignment.scan.min_cov,
    vircov_scan_reads: params.taxa.alignment.scan.min_reads,
    vircov_scan_coverage: params.taxa.alignment.scan.min_coverage,
    vircov_scan_regions: params.taxa.alignment.scan.min_regions,
    vircov_scan_regions_coverage: params.taxa.alignment.scan.min_mapq,
    vircov_group_by: params.taxa.alignment.scan.selection.group_by,
    vircov_group_sep: params.taxa.alignment.scan.selection.group_sep,
    vircov_group_select_by: params.taxa.alignment.scan.selection.select_by,
    vircov_group_select_segment_field: params.taxa.alignment.scan.selection.segment_field,
    vircov_group_select_segment_field_nan: params.taxa.alignment.scan.selection.segment_field_nan,
);
include { VircovRealign as VircovRealign } from '../../processes/vircov' addParams(
    subdir: "pathogens/alignment/remap",
    vircov_remap_min_mapq: params.taxa.alignment.remap.min_mapq,
    vircov_remap_min_len: params.taxa.alignment.remap.min_len,
    vircov_remap_min_cov: params.taxa.alignment.remap.min_cov,
    vircov_remap_reads: params.taxa.alignment.remap.min_reads,
    vircov_remap_coverage: params.taxa.alignment.remap.min_coverage,
    vircov_remap_regions: params.taxa.alignment.remap.min_regions,
    vircov_remap_regions_coverage: params.taxa.alignment.remap.min_mapq,
);
include { ConcatenateVircov as ConcatenateVircov } from '../../processes/vircov' addParams(
    subdir: "pathogens/alignment/remap"
);
include { ScrubbyExtractVircovReads as ScrubbyExtractVircovReads } from '../../processes/scrubby' addParams(
    subdir: "pathogens/alignment/remap"
);
include { MashScreenWinner as MashScreen } from '../../processes/mash' addParams(
    subdir: "pathogens/alignment/mash_subset"
);
include { MashDatabaseSubset as MashDatabaseSubset } from '../../processes/cerebro' addParams(
    subdir: "pathogens/alignment/mash_subset",
    min_shared_hashes: params.taxa.alignment.subset.min_shared_hashes
);
include { MinimapIndexSubset as MinimapIndexSubset } from '../../processes/minimap2' addParams(
    subdir: "pathogens/alignment/mash_subset"
);

/* 
==============================
DATABASE SUBSETS FOR ALIGNMENT
==============================
*/


workflow background {
    take: 
        reads
        references                                   
        kraken_dbs    
        domain                                             
    main: 
        output = ScrubbyBackgroundDepletion(reads, kraken_dbs, references)     
    emit: 
        reads = output.reads
        results = output.results                 
}

workflow subset {
    take:
        reads                                                                                     
        mash_index                                                                    
        db_fasta
    main:
        mash_screen = MashScreen(reads, mash_index)                              
        mash_subset = MashDatabaseSubset(mash_screen, db_fasta)       
        subset_db = MinimapIndexSubset(mash_subset)    
    emit:
        reads = reads
        index = subset_db.index
        fasta = subset_db.fasta
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
        alignment = MinimapAlign(reads, db_index)   
        reference_selection = VircovReferenceSelection(alignment, db_fasta, blacklist)
        // For each sample, emit the selected references with their input reads individually for parallelisation
        // If there is just one output from the dynamic reference file outputs, it needs to be put into an
        // array - this does not automatically happen - see the additional mapping step
        aligned_references = reference_selection.references | 
            map { data -> data[4] instanceof Collection ? data : tuple(data[0], data[1], data[2], data[3], [data[4]]) } |
            map { it[4].collect { el -> tuple(it[0], it[1], it[2], it[3], el) } } | 
            flatMap
    emit:
        references = aligned_references         // individual alignments and reference selections for paralellisation of realignments
        results = reference_selection.results
}

workflow remap {
    take: 
        aligned_references      
        domain           
    main: 
        // For each selected reference, realign the input reads and filter by alignment and coverage metrics
        realignment = MinimapRealign(aligned_references) 
        realignment_statistics = VircovRealign(realignment.reads, realignment.aligned)
        // Group outputs by sample and database name to concatenate results
        collected_statistics = realignment_statistics.aligned | groupTuple(by: [0, 1]) | ConcatenateVircov
        // ... and extract aligned reads for each reference 
        extracted_reads = realignment_statistics | ScrubbyExtractVircovReads
    emit:
        aligned = extracted_reads.aligned           // id, db_name, ref_name, ref_fasta, bam, vircov, fwd_extracted, rev_extracted
        coverage = realignment.coverage             // id, db_name, covtxt
        results = collected_statistics.results
}


workflow alignment_illumina {
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
        background(reads, references, kraken_dbs, domain)

        if (database_subset) {
            subset(background.out.reads, mash_index, db_fasta)
            scan(subset.out.reads, subset.out.index, subset.out.fasta, blacklist, domain)
            remap(scan_ont.out.references, domain)
        } else {
            scan(background.out.reads, db_index, db_fasta, blacklist, domain)
            remap(scan.out.references, domain)
        }
        
    emit:
        scan_results = scan.out.results
        remap_results = remap.out.results
        remap_aligned = remap.out.aligned
        remap_coverage = remap.out.coverage
}   
