

include { ScrubbyReadsKrakenMinimapDepletion as ScrubbyBackgroundDepletion } from '../../processes/scrubby' addParams(
    result_file: "qc__scrubby__background",
    subdir: "quality_control/background",
    kraken_taxa: params.qc.background.depletion.taxa,
    kraken_taxa_direct: params.qc.background.depletion.direct,
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
);

include { ScrubbyReadsKrakenMinimapDepletionOnt as ScrubbyBackgroundDepletionOnt } from '../../processes/scrubby' addParams(
    result_file: "qc__scrubby__background",
    subdir: "quality_control/background",
    kraken_taxa: params.qc.background.depletion.taxa,
    kraken_taxa_direct: params.qc.background.depletion.direct,
    deplete_min_cov: 0, 
    deplete_min_len: 0, 
    deplete_min_mapq: 0
);

workflow background_depletion {
    take: 
        reads
        references                                   
        kraken_dbs    
        ont                                                  
    main: 
        if (ont) {
            output = ScrubbyBackgroundDepletionOnt(reads, kraken_dbs, references)    
        } else {
            output = ScrubbyBackgroundDepletion(reads, kraken_dbs, references)    
        }
    emit: 
        reads = output.reads
        results = output.results                 
}