/* 
=========================
K-MER PROFILING WORKFLOWS
=========================
*/

include { Kraken2Uniq as Kraken2Uniq } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2uniq",
    minimum_hit_groups: params.kraken2_minimum_hit_groups
)

workflow kraken2uniq {
    take:
        reads      
        kraken_dbs      
    main:
        Kraken2Uniq(reads, kraken_dbs)
    emit:
        reads = reads
        results = Kraken2Uniq.out.results
}

// For appropriate level in run outputs
workflow kmer_profiling {
    take:
        reads
        kraken_dbs
    main:
        kraken2uniq(reads, kraken_dbs)
    emit:
        reads = reads
        results = kraken2uniq.out.results
}