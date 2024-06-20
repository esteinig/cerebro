/* 
=========================
K-MER PROFILING WORKFLOWS
=========================
*/

include { Kraken2Uniq } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2uniq",
    kraken2_minimum_hit_groups: params.taxa.kmer.kraken2.minimum_hit_groups
);
include { Kraken2Bracken } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2bracken",
    kraken2_confidence: params.taxa.kmer.kraken2.confidence,
    bracken_read_length: params.taxa.kmer.bracken.read_length,
    bracken_taxonomic_level: params.taxa.kmer.bracken.taxonomic_level,
    bracken_read_threshold: params.taxa.kmer.bracken.read_threshold
);
include { Metabuli } from '../../processes/metabuli' addParams(
    subdir: "kmer/metabuli"
);

include { Kraken2UniqOnt } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2uniq",
    kraken2_minimum_hit_groups: params.taxa.kmer.kraken2.minimum_hit_groups
);
include { Kraken2BrackenOnt } from '../../processes/kraken2' addParams(
    subdir: "kmer/kraken2bracken",
    kraken2_confidence: params.taxa.kmer.kraken2.confidence,
    bracken_read_length: params.taxa.kmer.bracken.read_length,
    bracken_taxonomic_level: params.taxa.kmer.bracken.taxonomic_level,
    bracken_read_threshold: params.taxa.kmer.bracken.read_threshold
);
include { MetabuliOnt } from '../../processes/metabuli' addParams(
    subdir: "kmer/metabuli"
)


workflow kraken2uniq {
    take:
        reads      
        kraken_dbs
        ont
    main:
        if (ont) {
            kraken2uniq = Kraken2UniqOnt(reads, kraken_dbs)
        } else {
            kraken2uniq = Kraken2Uniq(reads, kraken_dbs)
        }
        
    emit:
        reads = reads
        results = kraken2uniq[0]
}

workflow kraken2bracken {
    take:
        reads      
        kraken_dbs  
        ont    
    main:
        if (ont) {
            kraken2bracken = Kraken2BrackenOnt(reads, kraken_dbs)
        } else {
            kraken2bracken = Kraken2Bracken(reads, kraken_dbs)
        }
    emit:
        reads = reads
        results = kraken2bracken[0]
}


workflow metabuli {
    take:
        reads      
        metabuli_dbs  
        ont    
    main:
        if (ont) {
            metabuli = MetabuliOnt(reads, kraken_dbs)
        } else {
            metabuli = Metabuli(reads, kraken_dbs)
        }
    emit:
        reads = reads
        results = metabuli[0]
}

workflow kmer_profiling {
    take:
        reads
        kraken_dbs
        metabuli_dbs
        ont
    main:
        if (params.taxa.kmer.kraken2bracken.enabled) {
            output = kraken2bracken(reads, kraken_dbs, ont)
            kracken2bracken_results = output.results
        } else {
            kracken2bracken_results = Channel.empty()
        }
        if (params.taxa.kmer.metabuli.enabled) {
            output = metabuli(reads, metabuli_dbs, ont)
            metabuli_results = output.results
        } else {
            metabuli_results = Channel.empty()
        }

        results = kraken2bracken_results.mix(metabuli_results)

    emit:
        reads = reads
        results = results
}

