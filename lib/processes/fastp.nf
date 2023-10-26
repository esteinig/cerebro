process Fastp {

    label "fastp"
    tag { id }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "qc__fastp__reads" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"
    publishDir "$params.outdir/workflow/$params.subdir/fastq", mode: "symlink", pattern: "${id}_qc_*.fq.gz" 

    input:
    tuple val(id), path(forward), path(reverse)
    path(adapter_fasta)
    val(deduplicate)

    output:
    tuple (val(id), path("${id}_qc_1.fq.gz"), path("${id}_qc_2.fq.gz"), emit: reads)
    tuple (val(id), path("qc__fastp__reads"), emit: results)
    path("${id}.json")

    script:

    dedup = deduplicate ? "--dedup" : ""
    adapter_file = adapter_fasta ? "--adapter_fasta $adapter_fasta" : ""
    
    adapter_auto_detection = params.qc.reads.fastp.adapter_auto_detect ? "--detect_adapter_for_pe" : ""
    read_length = params.qc.reads.fastp.min_read_length ? "--length_required $params.qc.reads.fastp.min_read_length" : ""
    trim_poly_g = params.qc.reads.fastp.trim_poly_g ? " --trim_poly_g --poly_g_min_len $params.qc.reads.fastp.trim_poly_g" : ""
    tail_quality = params.qc.reads.fastp.cut_tail_quality ? "--cut_tail --cut_tail_mean_quality $params.qc.reads.fastp.cut_tail_quality" : ""
    complexity_filter = params.qc.reads.fastp.complexity_threshold ? "--low_complexity_filter --complexity_threshold $params.qc.reads.fastp.complexity_threshold" : ""   
    adapter_sequence = params.qc.reads.fastp.adapter_seq_1 && params.qc.reads.fastp.adapter_seq_2 ? "--adapter_sequence=$params.qc.reads.fastp.adapter_seq_1 --adapter_sequence_r2=$params.qc.reads.fastp.adapter_seq_2" : ""
    

    """
    fastp -i $forward -I $reverse -o ${id}_qc_1.fq.gz -O ${id}_qc_2.fq.gz --thread $task.cpus --json ${id}.json --dup_calc_accuracy 6  \
        $read_length $tail_quality $complexity_filter $adapter_auto_detection $adapter_file $adapter_sequence $dedup $trim_poly_g
    cp ${id}.json qc__fastp__reads
    """

}

// Scan reads only for total counts when using deduplication

process FastpScan {

    label "fastp"
    tag { id }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "qc__fastp__scan" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"

    input:
    tuple val(id), path(forward), path(reverse)

    output:
    tuple (val(id), path("qc__fastp__scan"), emit: results)
    path("${id}.json")

    script:

    """
    fastp -i $forward -I $reverse --thread $task.cpus --json ${id}.json \
        --disable_adapter_trimming --dont_eval_duplication --disable_trim_poly_g \
        --disable_quality_filtering --disable_length_filtering 
    cp ${id}.json qc__fastp__scan
    """

}