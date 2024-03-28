/* 
=========================
K-MER PROFILING WORKFLOWS
=========================
*/

include { Kraken2Uniq } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2uniq",
    kraken2_minimum_hit_groups: params.taxa.kmer.kraken2.minimum_hit_groups
)

include { Kraken2Bracken } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2bracken",
    kraken2_confidence: params.taxa.kmer.kraken2.confidence,
    bracken_read_length: params.taxa.kmer.bracken.read_length,
    bracken_taxonomic_level: params.taxa.kmer.bracken.taxonomic_level,
    bracken_read_threshold: params.taxa.kmer.bracken.read_threshold
)


include { Kraken2UniqOnt } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2uniq",
    kraken2_minimum_hit_groups: params.taxa.kmer.kraken2.minimum_hit_groups
)

include { Kraken2BrackenOnt } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2bracken",
    kraken2_confidence: params.taxa.kmer.kraken2.confidence,
    bracken_read_length: params.taxa.kmer.bracken.read_length,
    bracken_taxonomic_level: params.taxa.kmer.bracken.taxonomic_level,
    bracken_read_threshold: params.taxa.kmer.bracken.read_threshold
)



workflow kraken2uniq {
    take:
        reads      
        kraken_dbs
        ont
    main:
        if (ont) {
            k2u = Kraken2UniqOnt(reads, kraken_dbs)
        } else {
            k2u = Kraken2Uniq(reads, kraken_dbs)
        }
        
    emit:
        reads = reads
        results = k2u[0]
}

workflow kraken2bracken {
    take:
        reads      
        kraken_dbs  
        ont    
    main:
        if (ont) {
            k2b = Kraken2Bracken(reads, kraken_dbs)
        } else {
            k2b = Kraken2Bracken(reads, kraken_dbs)
        }
    emit:
        reads = reads
        results = k2b[0]
}


// Kraken2Uniq
workflow kmer_pathogen_detection {
    take:
        reads
        kraken_dbs
        ont
    main:
        kraken2uniq(reads, kraken_dbs, ont)
    emit:
        reads = reads
        results = kraken2uniq.out.results
}


// Kraken2Bracken
workflow kmer_pathogen_profiling {
    take:
        reads
        kraken_dbs
        ont
    main:
        kraken2bracken(reads, kraken_dbs, ont)
    emit:
        reads = reads
        results = kraken2bracken.out.results
}