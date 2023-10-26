/* 
====================
Taxonomic depletion
====================
*/

include { ScrubbyReadsKrakenMinimapDepletion as ScrubbyVirusBackgroundDepletion } from '../../processes/scrubby' addParams(
    subdir: "viruses/depletion/background",
    result_file: "qc__scrubby__virus_background",
    kraken_taxa: params.virus_background_taxa,
    kraken_taxa_direct: params.virus_background_direct,
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
)
include { MinimapAlignPAF as MinimapAlign } from '../../processes/minimap2' addParams(
    subdir: "viruses/scan/alignments"
)
include { MinimapRealignBAM as MinimapRealign } from '../../processes/minimap2' addParams(
    subdir: "viruses/remap/alignments",
    covtobed_max_cov: params.covtobed_max_cov
)
include { VircovReferenceSelection } from '../../processes/vircov' addParams(
    subdir: "viruses/scan",
    vircov_scan_min_mapq: params.vircov_scan_min_mapq,
    vircov_scan_min_len: params.vircov_scan_min_len,
    vircov_scan_min_cov: params.vircov_scan_min_cov,
    vircov_scan_reads: params.vircov_scan_reads,
    vircov_scan_coverage: params.vircov_scan_coverage,
    vircov_scan_regions: params.vircov_scan_regions,
    vircov_scan_regions_coverage: params.vircov_scan_regions_coverage,
    vircov_group_by: params.vircov_group_by,
    vircov_group_sep: params.vircov_group_sep,
    vircov_group_select_by: params.vircov_group_select_by

)
include { VircovRealign } from '../../processes/vircov' addParams(
    subdir: "viruses/remap",
    vircov_remap_min_mapq: params.vircov_remap_min_mapq,
    vircov_remap_min_len: params.vircov_remap_min_len,
    vircov_remap_min_cov: params.vircov_remap_min_cov,
    vircov_remap_reads: params.vircov_remap_reads,
    vircov_remap_coverage: params.vircov_remap_coverage,
    vircov_remap_regions: params.vircov_remap_regions,
    vircov_remap_regions_coverage: params.vircov_remap_regions_coverage
)
include { ConcatenateVircov } from '../../processes/vircov' addParams(
    subdir: "viruses/remap"
)
include { ScrubbyExtractVircovReads } from '../../processes/scrubby' addParams(
    subdir: "viruses/remap/reads"
)
include { IvarConsensus } from '../../processes/ivar' addParams(
    subdir: "viruses/assembly/consensus"
)
include { ProcessVirusAlignment } from '../../processes/cerebro' addParams(
    subdir: "viruses/"
)


workflow background_depletion {
    take: 
        reads
        references                                   
        kraken_dbs                                                       
    main: 
        ScrubbyVirusBackgroundDepletion(reads, kraken_dbs, references)       
    emit: 
        reads = ScrubbyVirusBackgroundDepletion.out.reads
        results = ScrubbyVirusBackgroundDepletion.out.results                 
}


workflow virus_scan {
    take: 
        reads                                                 
        db_index
        db_fasta
        blacklist                                    
    main: 
        // Align reads against viral database and select optimal references
        MinimapAlign(reads, db_index)   
        VircovReferenceSelection(MinimapAlign.out, db_fasta, blacklist)
        // For each sample, emit the selected references with their input reads individually for parallelisation
        // If there is just one output from the dynamic reference file outputs, it needs to be put into an
        // array - this does not automatically happen - see the additional mapping step
        aligned_references = VircovReferenceSelection.out.references |
            map { data -> data[4] instanceof Collection ? data : tuple(data[0], data[1], data[2], data[3], [data[4]]) } |
            map { it[4].collect { el -> tuple(it[0], it[1], it[2], it[3], el) } } | 
            flatMap
    emit:
        references = aligned_references // individual alignments and reference selections for paralellisation of realignments
        results = VircovReferenceSelection.out.results
}

workflow virus_remap {
    take: 
        aligned_references                                      
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
        coverage = MinimapRealign.out.coverage // id, db_name, covtxt
        results = ConcatenateVircov.out.results
}


workflow virus_assembly {
    take:
        aligned_data
    main:
        // Create a consensus reference sequence from aligned reads with iVar
        IvarConsensus(aligned_data)
    emit:
        IvarConsensus.out.consensus
}

workflow virus_detection {
    take:
        reads
        references
        kraken_dbs
        db_index
        db_fasta
        blacklist
        assembly
    main:

        // Background depletion
        background_depletion(reads, references, kraken_dbs)
        
        // Scan and remap 
        virus_scan(background_depletion.out.reads, db_index, db_fasta, blacklist) 
        virus_remap(virus_scan.out.references)

        // Optional consensus assembly
        if (assembly) { 
            virus_assembly(virus_remap.out.aligned)
            consensus = virus_assembly.out
        } else {
            consensus = Channel.empty()
        }

        // Scan-remap-consensus result aggregation
        result_tables = virus_scan.out.results.mix(virus_remap.out.results) | groupTuple(by: [0, 1])

        // Mixing outputs does not guarantee order, remap results are not guaranteed, 
        // we will account for this inside the process
        result_data = result_tables.mix(
            consensus | groupTuple(by: [0, 1]),
            virus_remap.out.coverage | groupTuple(by: [0, 1])
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
        ProcessVirusAlignment(result_data)

    emit:
        results = ProcessVirusAlignment.out.results
}   
