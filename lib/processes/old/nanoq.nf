
process Nanoq {

    label "nanoq"
    tag { id }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "qc__nanoq__reads" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "symlink", pattern: "${id}_qc.fq.gz"
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"

    input:
    tuple val(id), path(reads)

    output:
    tuple (val(id), path("${id}_qc.fq.gz"), emit: reads)
    tuple (val(id), path("qc__nanoq__reads"), emit: results)
    path("${id}.json")

    script:

    """
    nanoq -i $reads -l $params.qc.reads.nanoq.min_read_length -q $params.qc.reads.nanoq.min_read_quality --report ${id}.json --json > ${id}_qc.fq.gz
    cp ${id}.json qc__nanoq__reads
    """

}


process NanoqScan {

    label "nanoq"
    tag { id }

    publishDir "$params.outdir/results/$id", mode: "copy", pattern: "${id}.json", saveAs: { "qc__nanoq__scan" }
    publishDir "$params.outdir/workflow/$params.subdir", mode: "copy", pattern: "${id}.json"

    input:
    tuple val(id), path(reads)

    output:
    tuple (val(id), path("qc__nanoq__scan"), emit: results)
    path("${id}.json")

    script:

    """
    nanoq -i $reads --report ${id}.json --json -s
    cp ${id}.json qc__nanoq__scan
    """

}