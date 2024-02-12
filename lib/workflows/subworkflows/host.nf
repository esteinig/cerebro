
include { ScrubbyReadsKrakenMinimapDepletion as ScrubbyIllumina } from '../../processes/scrubby' addParams(
    subdir: "quality_control/host",
    result_file: "qc__scrubby__host",
    kraken_taxa: params.qc.host.depletion.taxa,
    kraken_taxa_direct: params.qc.host.depletion.direct,
    deplete_min_cov: params.qc.host.depletion.min_cov, 
    deplete_min_len: params.qc.host.depletion.min_len, 
    deplete_min_mapq: params.qc.host.depletion.min_mapq
)

include { ScrubbyReadsKrakenMinimapDepletionOnt as ScrubbyOnt } from '../../processes/scrubby' addParams(
    subdir: "quality_control/host",
    result_file: "qc__scrubby__host",
    kraken_taxa: params.qc.host.depletion.taxa,
    kraken_taxa_direct: params.qc.host.depletion.direct,
    deplete_min_cov: params.qc.host.depletion.min_cov, 
    deplete_min_len: params.qc.host.depletion.min_len, 
    deplete_min_mapq: params.qc.host.depletion.min_mapq
)

workflow host_depletion {
    take:
        reads
        kraken_dbs
        host_references
        ont
    main:
        if (ont) {
            scrubby = ScrubbyOnt(reads, kraken_dbs, host_references)
        } else {
            scrubby = ScrubbyIllumina(reads, kraken_dbs, host_references)
        }
    emit:
        reads = scrubby.reads
        host_results = scrubby.results 
}