
include { c } from '../../utils';

include { UmiToolsAlignment } from '../../processes/umi';
include { NaiveDeduplication } from '../../processes/umi';
include { CalibDeduplication } from '../../processes/umi';
include { UmiToolsDeduplication } from '../../processes/umi';

include { ercc_control } from './ercc';
include { ercc_control as ercc_clean } from './ercc';
include { phage_control } from './phage';

include { Fastp as ReadQualityControl } from '../../processes/fastp' addParams(
    subdir: "quality_control/read_qc",
)
include { FastpScan as ReadScan } from '../../processes/fastp' addParams(
    subdir: "quality_control/read_qc/scan",
)

include { ScrubbyReadsKrakenMinimapDepletion as HostDepletion } from '../../processes/scrubby' addParams(
    subdir: "quality_control/host",
    result_file: "qc__scrubby__host",
    kraken_taxa: params.qc.host.depletion.taxa,
    kraken_taxa_direct: params.qc.host.depletion.direct,
    deplete_min_cov: params.qc.host.depletion.min_cov, 
    deplete_min_len: params.qc.host.depletion.min_len, 
    deplete_min_mapq: params.qc.host.depletion.min_mapq
)

workflow quality_control {
    take:
        reads
        adapter_fasta
        ercc_fasta
        phage_fasta
        kraken_dbs
        host_references
        host_ercc_index
        deduplication
        deduplication_method
        host_removal
        phage_removal
    main:

        if (deduplication) {
            // If we deduplicate, we must scan the reads first 
            // to get the total read counts for processing later
            ReadScan(reads)
            scan_results = ReadScan.out.results

            if (deduplication_method == "umi-tools" || deduplication_method == "umi-tools-naive") {
                UmiToolsAlignment(reads, host_ercc_index, "host_ercc")
                reads = UmiToolsDeduplication(UmiToolsAlignment.out.alignment, params.qc.deduplication.seed)
                if (deduplication == "umi-tools-naive") {
                    reads = NaiveDeduplication(reads, true)      
                }
            } else if (deduplication_method == "umi-calib"){
                reads = CalibDeduplication(reads)
            } else if (deduplication_method == "naive"){
                reads = NaiveDeduplication(reads, false)
            } else if (deduplication_method == "umi-naive"){
                reads = NaiveDeduplication(reads, true)
            }
        } else {
            scan_results = Channel.empty()
        }
        
        // ERCCs are aligned and removed
        if (params.qc.controls.ercc.enabled) {

            ercc_control(reads, ercc_fasta)
            reads = ercc_control.out.reads
            ercc_results = ercc_control.out.results

        } else {
            ercc_results = Channel.empty()
        }

        // Read quality control and removal of low-complexity reads,
        // optional de-duplication of identical sequences
        ReadQualityControl(reads, adapter_fasta, false)
        reads = ReadQualityControl.out.reads
        qc_results = ReadQualityControl.out.results

        // Sequential k-mer and alignment host depletion
        if (host_removal) {

            HostDepletion(reads, kraken_dbs, host_references)
            reads = HostDepletion.out.reads
            host_results = HostDepletion.out.results

        } else {
            host_results = Channel.empty()
        }

        // Phage spike-ins are aligned and removed
        // after human removal for faster read alignments
        if (phage_removal) {

            phage_control(reads, phage_fasta)
            reads = phage_control.out.reads
            phage_results = phage_control.out.results

        } else {
            phage_results = Channel.empty()
        }

    emit:
        reads = reads
        results = qc_results.mix(scan_results, ercc_results, host_results, phage_results)
}


