
include { c } from '../../utils';

include { ercc_control } from './ercc';
include { phage_control } from './phage';
include { host_depletion } from './host';
include { background_depletion } from './background';

include { NaiveDeduplication } from '../../processes/umi';
include { CalibDeduplication } from '../../processes/umi';

include { Nanoq } from '../../processes/nanoq' addParams(
    subdir: "quality_control/read_qc",
);
include { NanoqScan } from '../../processes/nanoq' addParams(
    subdir: "quality_control/read_qc/scan",
);
include { Fastp } from '../../processes/fastp' addParams(
    subdir: "quality_control/read_qc",
);
include { FastpScan } from '../../processes/fastp' addParams(
    subdir: "quality_control/read_qc/scan",
);

workflow quality_control_illumina {
    take:
        reads // id, r1, r2
        adapter_fasta
        ercc_fasta
        phage_fasta
        host_kraken_dbs
        host_references
        background_kraken_dbs
        background_references
        host_removal
        phage_removal
        background_removal
        deduplication
        deduplication_method
    main:
        // If we deduplicate, we must scan the reads 
        // to get the total read counts for parsing
        FastpScan(reads)
        scan_results = FastpScan.out.results

        if (deduplication) {
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
        // optional de-duplication of identical sequences (Fastp) but
        // note that unlike native de-duplication the algorithm in
        // Fastp has a relatively high chance of hash-collision, so
        // is essentially non-deterministic

        if (params.qc.reads.fastp.enabled) {
            Fastp(reads, adapter_fasta, false)
            reads = Fastp.out.reads
            qc_results = Fastp.out.results
        } else {
            qc_results = Channel.empty()
        }
        

        // Sequential k-mer and alignment host depletion
        if (host_removal) {
            host_depletion(reads, host_kraken_dbs, host_references, false)
            reads = host_depletion.out.reads
            host_results = host_depletion.out.results
        } else {
            host_results = Channel.empty()
        }

        // Phage spike-ins are aligned and removed
        if (phage_removal) {
            phage_control(reads, phage_fasta, false)
            reads = phage_control.out.reads
            phage_results = phage_control.out.results
        } else {
            phage_results = Channel.empty()
        }

        // Other background removal is conducted
        if (background_removal) {
            background_depletion(reads, background_kraken_dbs, background_references, false)
            reads = background_depletion.out.reads
            background_results = background_depletion.out.results
        } else {
            background_results = Channel.empty()
        }

    emit:
        reads = reads
        results = qc_results.mix(scan_results, ercc_results, host_results, phage_results, background_results)
}



workflow quality_control_ont {
    take:
        reads  // id, fq
        ercc_fasta
        phage_fasta
        host_kraken_dbs
        host_references
        background_kraken_dbs
        background_references
        host_removal
        phage_removal
        background_removal
    main:
    
        NanoqScan(reads)
        scan_results = NanoqScan.out.results

        if (params.qc.controls.ercc.enabled) {
            ercc = ercc_control(reads, ercc_fasta, true)
            reads = ercc.reads
            ercc_results = ercc.results
        } else {
            ercc_results = Channel.empty()
        }
        
        if (params.qc.reads.fastp.enabled) {
            nanoq = Nanoq(reads)
            reads = Nanoq.out.reads
            qc_results = Nanoq.out.results
        } else {
            qc_results = Channel.empty()
        }
        
        // Sequential k-mer and alignment host depletion
        if (host_removal) {
            host = host_depletion(nanoq.reads, host_kraken_dbs, host_references, true)
            reads = host.reads
            host_results = host.results
        } else {
            host_results = Channel.empty()
        }

        // Phage spike-ins are aligned and removed
        if (phage_removal) {
            phage_control = phage_control(reads, phage_fasta, true)

            reads = phage_control.reads
            phage_results = phage_control.results
        } else {
            phage_results = Channel.empty()
        }

        // Other background removal
        if (background_removal) {
            background_depletion(reads, background_kraken_dbs, background_references, true)
            reads = background_depletion.out.reads
            background_results = background_depletion.out.results
        } else {
            background_results = Channel.empty()
        }

    emit:
        reads = reads
        results = qc_results.mix(scan_results, ercc_results, host_results, phage_results, background_results)
}

