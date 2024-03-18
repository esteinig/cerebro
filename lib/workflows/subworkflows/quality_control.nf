
include { c } from '../../utils';

include { NaiveDeduplication } from '../../processes/umi';
include { CalibDeduplication } from '../../processes/umi';

include { ercc_control } from './ercc';
include { phage_control } from './phage';


include { Nanoq } from '../../processes/nanoq' addParams(
    subdir: "quality_control/read_qc",
)
include { Fastp } from '../../processes/fastp' addParams(
    subdir: "quality_control/read_qc",
)
include { FastpScan } from '../../processes/fastp' addParams(
    subdir: "quality_control/read_qc/scan",
)

workflow quality_control_illumina {
    take:
        reads // id, r1, r2
        adapter_fasta
        ercc_fasta
        phage_fasta
        kraken_dbs
        host_references
        deduplication
        deduplication_method
        host_removal
        phage_removal
    main:

        if (deduplication) {
            // If we deduplicate, we must scan the reads 
            // to get the total read counts for parsing
            FastpScan(reads)
            scan_results = FastpScan.out.results

            if (deduplication_method == "umi-calib"){
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
            ercc_control(reads, ercc_fasta, false)
            
            reads = ercc_control.out.reads
            ercc_results = ercc_control.out.results
        } else {
            ercc_results = Channel.empty()
        }

        // Read quality control and removal of low-complexity reads,
        // optional de-duplication of identical sequences
        Fastp(reads, adapter_fasta, false)
        reads = Fastp.out.reads
        qc_results = Fastp.out.results

        // Sequential k-mer and alignment host depletion
        if (host_removal) {
            host_depletion(reads, kraken_dbs, host_references, false)
            reads = host_depletion.out.reads
            host_results = host_depletion.out.results
        } else {
            host_results = Channel.empty()
        }

        // Phage spike-ins are aligned and removed
        // after human removal for faster read alignments
        if (phage_removal) {
            phage_control(reads, phage_fasta, false)
            reads = phage_control.out.reads
            phage_results = phage_control.out.results

        } else {
            phage_results = Channel.empty()
        }

    emit:
        reads = reads
        results = qc_results.mix(scan_results, ercc_results, host_results, phage_results)
}



workflow quality_control_ont {
    take:
        reads  // id, fq
        ercc_fasta
        phage_fasta
        kraken_dbs
        host_references
        host_removal
        phage_removal
    main:
        
        if (params.qc.controls.ercc.enabled) {
            ercc_control(reads, ercc_fasta, true)
            reads = ercc_control.out.reads
            ercc_results = ercc_control.out.results
        } else {
            ercc_results = Channel.empty()
        }

        Nanoq(reads)
        reads = Nanoq.out.reads
        qc_results = Nanoq.out.results

        if (host_removal) {
            host_depletion(reads, kraken_dbs, host_references, true)
            reads = host_depletion.out.reads
            host_results = host_depletion.out.results
        } else {
            host_results = Channel.empty()
        }

        if (phage_removal) {
            phage_control(reads, phage_fasta, true)
            reads = phage_control.out.reads
            phage_results = phage_control.out.results
        } else {
            phage_results = Channel.empty()
        }

    emit:
        reads = reads
        results = qc_results.mix(ercc_results, host_results, phage_results)
}

