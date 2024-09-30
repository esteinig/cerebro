include { MinimapAlignPaf as MinimapAlign } from '../../processes/minimap2' addParams(
    subdir: "pathogens/alignment/scan"
);
include { MinimapRealignBam as MinimapRealign } from '../../processes/minimap2' addParams(
    subdir: "pathogens/alignment/remap"
);
include { VircovReferenceSelection as VircovReferenceSelection } from '../../processes/vircov' addParams(
    subdir: "pathogens/alignment/scan"
);
include { VircovRealign as VircovRealign } from '../../processes/vircov' addParams(
    subdir: "pathogens/alignment/remap"
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
        db_index
        db_fasta
        blacklist
        consensus_assembly
        database_subset
        mash_index
    main:

        if (database_subset) {
            subset(reads, mash_index, db_fasta)
            scan(subset.out.reads, subset.out.index, subset.out.fasta, blacklist)
            remap(scan_ont.out.references)
        } else {
            scan(reads, db_index, db_fasta, blacklist)
            remap(scan.out.references)
        }
        
    emit:
        scan_results = scan.out.results
        remap_results = remap.out.results
        remap_aligned = remap.out.aligned
        remap_coverage = remap.out.coverage
}   
