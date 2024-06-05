include { alignment_illumina } from '../subworkflows/alignment_illumina';
include { alignment_ont } from '../subworkflows/alignment_ont';

include { IvarConsensus } from '../../processes/ivar' addParams(
    subdir: "pathogens/alignment/consensus_assembly"
);
include { AlignmentTable } from '../../processes/cerebro' addParams(
    subdir: "pathogens/alignment"
);


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
        ont
    main:

        if (ont) {
            alignment = alignment_ont(
                reads,
                references,
                kraken_dbs,
                db_index,
                db_fasta,
                blacklist,
                consensus_assembly,
                database_subset,
                mash_index,
                domain
            )
        } else {
            alignment = alignment_illumina(
                reads,
                references,
                kraken_dbs,
                db_index,
                db_fasta,
                blacklist,
                consensus_assembly,
                database_subset,
                mash_index,
                domain
            )
        }
        
        // Optional consensus assembly
        if (consensus_assembly) { 
            consensus = virus_consensus_assembly(alignment.remap_aligned)
        } else {
            consensus = Channel.empty()
        }

        // Scan-remap-consensus result aggregation
        result_tables = alignment.scan_results.mix(alignment.remap_results) | groupTuple(by: [0, 1])

        // Mixing outputs does not guarantee order, remap results are not guaranteed, 
        // we will account for this inside the process
        result_data = result_tables.mix(
            consensus | groupTuple(by: [0, 1]),
            alignment.remap_coverage | groupTuple(by: [0, 1])
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
        AlignmentTable(result_data)
    emit:
        results = AlignmentTable.out.results
}   
